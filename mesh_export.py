#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""3D mesh export for molecule-icon-generator.

Builds triangulated meshes (UV sphere for atoms, capped cylinder for bonds)
from an RDKit molecule and writes them as ASCII FBX 7.4 — readable by Blender,
Maya, 3ds Max, Cinema 4D, Unity and Unreal.

Also includes a Wavefront OBJ + MTL writer as a dependency-free fallback.
"""

from __future__ import annotations

import math
import os
import time
from typing import Any, Iterable, List, Sequence, Tuple

import numpy as np
from scipy.spatial.transform import Rotation as Rot

# RDKit is only needed when consumers call build_molecule_meshes / save_3d_as_fbx.
# Importing it lazily keeps the mesh primitives and FBX writer usable on their
# own (e.g. for tests, or if a caller already has its own atom positions).


# ---------------------------------------------------------------------------
# Color helpers
# ---------------------------------------------------------------------------

def hex_to_rgb01(hex_color: str) -> Tuple[float, float, float]:
    """Convert '#rrggbb' to (r, g, b) floats in [0, 1]."""
    if hex_color is None:
        return 0.5, 0.5, 0.5
    h = hex_color.lstrip('#').strip()
    if len(h) != 6:
        return 0.5, 0.5, 0.5
    try:
        return (int(h[0:2], 16) / 255.0,
                int(h[2:4], 16) / 255.0,
                int(h[4:6], 16) / 255.0)
    except ValueError:
        return 0.5, 0.5, 0.5


# ---------------------------------------------------------------------------
# Primitive mesh builders
# ---------------------------------------------------------------------------

def sphere_mesh(center: Sequence[float], radius: float,
                segments: int = 24, rings: int | None = None
                ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Triangulated UV-sphere with smooth per-vertex normals.

    Returns (vertices [N,3], normals [N,3], triangles [M,3]).
    """
    if rings is None:
        rings = max(2, segments // 2)
    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
    cols = segments + 1  # closing column duplicates the seam
    n_verts = (rings + 1) * cols
    verts = np.empty((n_verts, 3), dtype=np.float32)
    norms = np.empty((n_verts, 3), dtype=np.float32)
    idx = 0
    for i in range(rings + 1):
        v = math.pi * i / rings
        sv = math.sin(v); cv = math.cos(v)
        for j in range(cols):
            u = 2.0 * math.pi * j / segments
            su = math.sin(u); cu = math.cos(u)
            nx = sv * cu; ny = sv * su; nz = cv
            verts[idx] = (cx + radius * nx, cy + radius * ny, cz + radius * nz)
            norms[idx] = (nx, ny, nz)
            idx += 1
    tris = np.empty((rings * segments * 2, 3), dtype=np.int32)
    t = 0
    for i in range(rings):
        for j in range(segments):
            a = i * cols + j
            b = a + 1
            c = a + cols
            d = c + 1
            tris[t] = (a, b, d); t += 1
            tris[t] = (a, d, c); t += 1
    return verts, norms, tris


def cylinder_mesh(start: Sequence[float], end: Sequence[float], radius: float,
                  segments: int = 24, capped: bool = True
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Triangulated cylinder between two points with optional flat caps.

    Side-vertex normals are radial (smooth side shading), cap-vertex normals
    point along the axis (flat caps). Returns (vertices, normals, triangles).
    """
    s = np.asarray(start, dtype=np.float64)
    e = np.asarray(end, dtype=np.float64)
    axis = e - s
    h = float(np.linalg.norm(axis))
    if h < 1e-9:
        empty = np.zeros((0, 3), dtype=np.float32)
        return empty, empty.copy(), np.zeros((0, 3), dtype=np.int32)
    az = axis / h
    # build orthonormal frame perpendicular to axis
    ref = np.array([1.0, 0.0, 0.0]) if abs(az[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(az, ref)
    u /= np.linalg.norm(u)
    v = np.cross(az, u)

    cols = segments + 1
    cos_t = np.cos(np.linspace(0.0, 2.0 * math.pi, cols))
    sin_t = np.sin(np.linspace(0.0, 2.0 * math.pi, cols))
    radial = (u[None, :] * cos_t[:, None] + v[None, :] * sin_t[:, None])  # (cols, 3)

    # side rings
    side_start = s[None, :] + radial * radius  # (cols, 3)
    side_end = e[None, :] + radial * radius

    side_verts = np.concatenate([side_start, side_end], axis=0).astype(np.float32)
    side_norms = np.concatenate([radial, radial], axis=0).astype(np.float32)

    side_tris = []
    for j in range(segments):
        a = j           # start ring j
        b = j + 1       # start ring j+1
        c = cols + j    # end ring j
        d = cols + j + 1
        side_tris.append((a, c, d))
        side_tris.append((a, d, b))

    if not capped:
        return side_verts, side_norms, np.array(side_tris, dtype=np.int32)

    # caps: each cap = 1 center vertex + segments+1 ring vertices, all sharing
    # the cap normal (-axis or +axis) so the caps shade flat.
    cap_norm_neg = (-az).astype(np.float32)
    cap_norm_pos = az.astype(np.float32)

    cap_a_verts = np.concatenate([s[None, :], side_start], axis=0).astype(np.float32)
    cap_a_norms = np.tile(cap_norm_neg, (cap_a_verts.shape[0], 1))
    cap_b_verts = np.concatenate([e[None, :], side_end], axis=0).astype(np.float32)
    cap_b_norms = np.tile(cap_norm_pos, (cap_b_verts.shape[0], 1))

    n_side = side_verts.shape[0]
    n_capa = cap_a_verts.shape[0]

    verts = np.concatenate([side_verts, cap_a_verts, cap_b_verts], axis=0)
    norms = np.concatenate([side_norms, cap_a_norms, cap_b_norms], axis=0)

    tris = list(side_tris)
    base_a = n_side          # cap A center index
    for j in range(segments):
        center_i = base_a
        ring_j = base_a + 1 + j
        ring_jp1 = base_a + 1 + j + 1
        # cap A normal points -axis, so wind clockwise as seen from outside
        tris.append((center_i, ring_jp1, ring_j))
    base_b = n_side + n_capa
    for j in range(segments):
        center_i = base_b
        ring_j = base_b + 1 + j
        ring_jp1 = base_b + 1 + j + 1
        tris.append((center_i, ring_j, ring_jp1))

    return verts, norms, np.array(tris, dtype=np.int32)


# ---------------------------------------------------------------------------
# Mesh container
# ---------------------------------------------------------------------------

class MeshPart:
    __slots__ = ('name', 'vertices', 'normals', 'triangles', 'color_hex')

    def __init__(self, name: str, vertices: np.ndarray, normals: np.ndarray,
                 triangles: np.ndarray, color_hex: str):
        self.name = name
        self.vertices = np.ascontiguousarray(vertices, dtype=np.float64)
        self.normals = np.ascontiguousarray(normals, dtype=np.float64)
        self.triangles = np.ascontiguousarray(triangles, dtype=np.int64)
        self.color_hex = color_hex


# ---------------------------------------------------------------------------
# Build meshes from an RDKit molecule
# ---------------------------------------------------------------------------

# Imported lazily inside the function to avoid a circular import.
def build_molecule_meshes(mol: Any,
                          atom_color: dict,
                          radius_multi: dict,
                          atom_radius: float = 0.5,
                          pos_multi: float = 1.0,
                          resolution: int = 24,
                          remove_H: bool = True,
                          rotation: Tuple[float, float, float] = (0, 0, 0),
                          bond_color: str | None = None
                          ) -> List[MeshPart]:
    """Walk the RDKit molecule and build a MeshPart per atom and per bond.

    Mirrors the geometry of `molecule_icon_generator.graph_3d` so the FBX
    export looks the same as the interactive plotly view.
    """
    from rdkit import Chem  # local import: only required for this entry point
    if remove_H:
        mol = Chem.RemoveHs(mol)
    if mol.GetNumConformers() == 0:
        raise ValueError("Molecule has no 3D conformer. Generate one first "
                         "(parse_structure(..., dimension_3=True)).")
    conf = mol.GetConformer()
    rotate = Rot.from_euler('xyz', rotation, degrees=True)

    # gather rotated, scaled positions
    pos = {}
    for i in range(mol.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pos[i] = rotate.apply((p.x, p.y, p.z)) * pos_multi

    bond_radius = atom_radius * radius_multi.get('Bond', 1.0) / 4.0
    bond_hex = bond_color if bond_color is not None else atom_color.get('Bond', '#444444')

    parts: List[MeshPart] = []

    # atoms
    for i in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        symbol = atom.GetSymbol()
        r = atom_radius * radius_multi.get(symbol, 1.0)
        color = atom_color.get(symbol, '#888888')
        v, n, t = sphere_mesh(pos[i], r, segments=max(8, resolution))
        parts.append(MeshPart(f'atom_{i}_{symbol}', v, n, t, color))

    # bonds
    for bond in mol.GetBonds():
        i1 = bond.GetBeginAtomIdx()
        i2 = bond.GetEndAtomIdx()
        v, n, t = cylinder_mesh(pos[i1], pos[i2], bond_radius,
                                segments=max(8, resolution))
        if v.shape[0] == 0:
            continue
        parts.append(MeshPart(f'bond_{bond.GetIdx()}_{i1}_{i2}', v, n, t,
                              bond_hex))

    return parts


# ---------------------------------------------------------------------------
# OBJ writer (handy fallback / sanity check)
# ---------------------------------------------------------------------------

def write_obj(parts: Sequence[MeshPart], filepath: str) -> str:
    """Write meshes as Wavefront OBJ + sibling MTL file."""
    base = os.path.splitext(filepath)[0]
    obj_path = base + '.obj'
    mtl_path = base + '.mtl'

    # one material per unique color
    unique_colors = []
    color_to_mat = {}
    for p in parts:
        if p.color_hex not in color_to_mat:
            mat_name = f'mat_{len(unique_colors)}'
            color_to_mat[p.color_hex] = mat_name
            unique_colors.append((mat_name, p.color_hex))

    with open(mtl_path, 'w') as mf:
        for mat_name, hex_color in unique_colors:
            r, g, b = hex_to_rgb01(hex_color)
            mf.write(f'newmtl {mat_name}\n')
            mf.write(f'Ka 0.0 0.0 0.0\n')
            mf.write(f'Kd {r:.4f} {g:.4f} {b:.4f}\n')
            mf.write(f'Ks 0.2 0.2 0.2\n')
            mf.write(f'Ns 20.0\n')
            mf.write(f'illum 2\n\n')

    with open(obj_path, 'w') as f:
        f.write(f'# molecule-icon-generator OBJ export\n')
        f.write(f'mtllib {os.path.basename(mtl_path)}\n')
        v_offset = 0
        n_offset = 0
        for p in parts:
            f.write(f'o {p.name}\n')
            f.write(f'usemtl {color_to_mat[p.color_hex]}\n')
            for v in p.vertices:
                f.write(f'v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n')
            for n in p.normals:
                f.write(f'vn {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n')
            f.write('s 1\n')
            for tri in p.triangles:
                a = int(tri[0]) + 1 + v_offset
                b = int(tri[1]) + 1 + v_offset
                c = int(tri[2]) + 1 + v_offset
                an = int(tri[0]) + 1 + n_offset
                bn = int(tri[1]) + 1 + n_offset
                cn = int(tri[2]) + 1 + n_offset
                f.write(f'f {a}//{an} {b}//{bn} {c}//{cn}\n')
            v_offset += p.vertices.shape[0]
            n_offset += p.normals.shape[0]
    return obj_path


# ---------------------------------------------------------------------------
# ASCII FBX 7.4 writer
# ---------------------------------------------------------------------------

def _safe(name: str) -> str:
    return ''.join(c if (c.isalnum() or c in '_-') else '_' for c in name)


_ID_COUNTER = [10_000_000]


def _fresh_id() -> int:
    _ID_COUNTER[0] += 17
    return _ID_COUNTER[0]


def _format_array(values: Iterable[float], per_line: int = 12,
                  precision: int = 6, indent: str = '\t\t\t\t') -> str:
    pieces = []
    n = 0
    fmt = f'{{:.{precision}f}}'
    for i, val in enumerate(values):
        sep = ',' if i else ''
        if n == per_line:
            pieces.append('\n' + indent)
            n = 0
        pieces.append(sep + fmt.format(float(val)))
        n += 1
    return ''.join(pieces)


def _format_int_array(values: Iterable[int], per_line: int = 24,
                      indent: str = '\t\t\t\t') -> str:
    pieces = []
    n = 0
    for i, val in enumerate(values):
        sep = ',' if i else ''
        if n == per_line:
            pieces.append('\n' + indent)
            n = 0
        pieces.append(sep + str(int(val)))
        n += 1
    return ''.join(pieces)


def write_ascii_fbx(parts: Sequence[MeshPart], filepath: str,
                    creator: str = 'molecule-icon-generator',
                    up_axis: str = 'Z') -> str:
    """Write meshes as an ASCII FBX 7.4 file at *filepath*.

    Materials are deduplicated by hex color so all carbons share one material,
    all hydrogens share another, etc.
    """
    if not filepath.lower().endswith('.fbx'):
        filepath = filepath + '.fbx'
    if not parts:
        raise ValueError('Cannot write FBX: no mesh parts supplied.')

    # one material per unique color
    color_to_mat: dict[str, dict] = {}
    for p in parts:
        if p.color_hex not in color_to_mat:
            color_to_mat[p.color_hex] = {
                'id': _fresh_id(),
                'name': f'mat_{len(color_to_mat)}',
                'color': p.color_hex,
            }

    geom_count = len(parts)
    model_count = len(parts)
    mat_count = len(color_to_mat)

    # assign ids — keep them in a side dict so MeshPart can use __slots__
    part_ids = {id(p): {'gid': _fresh_id(), 'mid': _fresh_id()} for p in parts}

    now = time.localtime()
    up_axis = up_axis.upper()
    up_index = {'X': 0, 'Y': 1, 'Z': 2}.get(up_axis, 2)

    out: List[str] = []
    a = out.append

    a('; FBX 7.4.0 project file')
    a(f'; Created by {creator}')
    a('; ----------------------------------------------------')
    a('')
    a('FBXHeaderExtension:  {')
    a('\tFBXHeaderVersion: 1003')
    a('\tFBXVersion: 7400')
    a('\tCreationTimeStamp:  {')
    a('\t\tVersion: 1000')
    a(f'\t\tYear: {now.tm_year}')
    a(f'\t\tMonth: {now.tm_mon}')
    a(f'\t\tDay: {now.tm_mday}')
    a(f'\t\tHour: {now.tm_hour}')
    a(f'\t\tMinute: {now.tm_min}')
    a(f'\t\tSecond: {now.tm_sec}')
    a('\t\tMillisecond: 0')
    a('\t}')
    a(f'\tCreator: "{creator}"')
    a('\tSceneInfo: "SceneInfo::GlobalInfo", "UserData" {')
    a('\t\tType: "UserData"')
    a('\t\tVersion: 100')
    a('\t\tMetaData:  {')
    a('\t\t\tVersion: 100')
    a('\t\t\tTitle: ""')
    a('\t\t\tSubject: ""')
    a('\t\t\tAuthor: ""')
    a('\t\t\tKeywords: ""')
    a('\t\t\tRevision: ""')
    a('\t\t\tComment: ""')
    a('\t\t}')
    a('\t}')
    a('}')
    a(f'Creator: "{creator}"')
    a('GlobalSettings:  {')
    a('\tVersion: 1000')
    a('\tProperties70:  {')
    a(f'\t\tP: "UpAxis", "int", "Integer", "",{up_index}')
    a('\t\tP: "UpAxisSign", "int", "Integer", "",1')
    a(f'\t\tP: "FrontAxis", "int", "Integer", "",{1 if up_index != 1 else 2}')
    a('\t\tP: "FrontAxisSign", "int", "Integer", "",1')
    a('\t\tP: "CoordAxis", "int", "Integer", "",0')
    a('\t\tP: "CoordAxisSign", "int", "Integer", "",1')
    a(f'\t\tP: "OriginalUpAxis", "int", "Integer", "",{up_index}')
    a('\t\tP: "OriginalUpAxisSign", "int", "Integer", "",1')
    a('\t\tP: "UnitScaleFactor", "double", "Number", "",1')
    a('\t\tP: "OriginalUnitScaleFactor", "double", "Number", "",1')
    a('\t\tP: "AmbientColor", "ColorRGB", "Color", "",0,0,0')
    a('\t\tP: "DefaultCamera", "KString", "", "", "Producer Perspective"')
    a('\t\tP: "TimeMode", "enum", "", "",11')
    a('\t}')
    a('}')
    a('Documents:  {')
    a('\tCount: 1')
    doc_id = _fresh_id()
    a(f'\tDocument: {doc_id}, "", "Scene" {{')
    a('\t\tProperties70:  {')
    a('\t\t\tP: "SourceObject", "object", "", ""')
    a('\t\t\tP: "ActiveAnimStackName", "KString", "", "", ""')
    a('\t\t}')
    a('\t\tRootNode: 0')
    a('\t}')
    a('}')
    a('References:  {')
    a('}')
    a('Definitions:  {')
    a('\tVersion: 100')
    a(f'\tCount: {1 + geom_count + model_count + mat_count}')
    a('\tObjectType: "GlobalSettings" {')
    a('\t\tCount: 1')
    a('\t}')
    a('\tObjectType: "Geometry" {')
    a(f'\t\tCount: {geom_count}')
    a('\t\tPropertyTemplate: "FbxMesh" {')
    a('\t\t\tProperties70:  {')
    a('\t\t\t\tP: "Color", "ColorRGB", "Color", "",0.8,0.8,0.8')
    a('\t\t\t\tP: "Primary Visibility", "bool", "", "",1')
    a('\t\t\t}')
    a('\t\t}')
    a('\t}')
    a('\tObjectType: "Model" {')
    a(f'\t\tCount: {model_count}')
    a('\t\tPropertyTemplate: "FbxNode" {')
    a('\t\t\tProperties70:  {')
    a('\t\t\t\tP: "Visibility", "Visibility", "", "A",1')
    a('\t\t\t}')
    a('\t\t}')
    a('\t}')
    a('\tObjectType: "Material" {')
    a(f'\t\tCount: {mat_count}')
    a('\t\tPropertyTemplate: "FbxSurfacePhong" {')
    a('\t\t\tProperties70:  {')
    a('\t\t\t\tP: "DiffuseColor", "ColorRGB", "Color", "",0.8,0.8,0.8')
    a('\t\t\t\tP: "DiffuseFactor", "Number", "", "A",1')
    a('\t\t\t\tP: "SpecularColor", "ColorRGB", "Color", "",0.2,0.2,0.2')
    a('\t\t\t\tP: "SpecularFactor", "Number", "", "A",0.5')
    a('\t\t\t\tP: "ShininessExponent", "Number", "", "A",20')
    a('\t\t\t}')
    a('\t\t}')
    a('\t}')
    a('}')
    a('Objects:  {')

    # Geometries
    for p in parts:
        gid = part_ids[id(p)]['gid']
        name = _safe(p.name)
        flat_verts = p.vertices.reshape(-1)
        flat_norms = p.normals.reshape(-1)
        # PolygonVertexIndex: last index of each polygon must be ~i (i.e., -i-1)
        idx = p.triangles.copy()
        idx[:, 2] = -idx[:, 2] - 1
        flat_idx = idx.reshape(-1)

        a(f'\tGeometry: {gid}, "Geometry::{name}", "Mesh" {{')
        a(f'\t\tVertices: *{flat_verts.size} {{')
        a('\t\t\ta: ' + _format_array(flat_verts))
        a('\t\t}')
        a(f'\t\tPolygonVertexIndex: *{flat_idx.size} {{')
        a('\t\t\ta: ' + _format_int_array(flat_idx))
        a('\t\t}')
        a('\t\tGeometryVersion: 124')
        a('\t\tLayerElementNormal: 0 {')
        a('\t\t\tVersion: 102')
        a('\t\t\tName: ""')
        a('\t\t\tMappingInformationType: "ByVertice"')
        a('\t\t\tReferenceInformationType: "Direct"')
        a(f'\t\t\tNormals: *{flat_norms.size} {{')
        a('\t\t\t\ta: ' + _format_array(flat_norms))
        a('\t\t\t}')
        a('\t\t}')
        a('\t\tLayerElementMaterial: 0 {')
        a('\t\t\tVersion: 101')
        a('\t\t\tName: ""')
        a('\t\t\tMappingInformationType: "AllSame"')
        a('\t\t\tReferenceInformationType: "IndexToDirect"')
        a('\t\t\tMaterials: *1 {')
        a('\t\t\t\ta: 0')
        a('\t\t\t}')
        a('\t\t}')
        a('\t\tLayer: 0 {')
        a('\t\t\tVersion: 100')
        a('\t\t\tLayerElement:  {')
        a('\t\t\t\tType: "LayerElementNormal"')
        a('\t\t\t\tTypedIndex: 0')
        a('\t\t\t}')
        a('\t\t\tLayerElement:  {')
        a('\t\t\t\tType: "LayerElementMaterial"')
        a('\t\t\t\tTypedIndex: 0')
        a('\t\t\t}')
        a('\t\t}')
        a('\t}')

    # Models (one per part)
    for p in parts:
        mid = part_ids[id(p)]['mid']
        name = _safe(p.name)
        a(f'\tModel: {mid}, "Model::{name}", "Mesh" {{')
        a('\t\tVersion: 232')
        a('\t\tProperties70:  {')
        a('\t\t\tP: "Visibility", "Visibility", "", "A",1')
        a('\t\t\tP: "Lcl Translation", "Lcl Translation", "", "A",0,0,0')
        a('\t\t\tP: "Lcl Rotation", "Lcl Rotation", "", "A",0,0,0')
        a('\t\t\tP: "Lcl Scaling", "Lcl Scaling", "", "A",1,1,1')
        a('\t\t\tP: "DefaultAttributeIndex", "int", "Integer", "",0')
        a('\t\t\tP: "InheritType", "enum", "", "",1')
        a('\t\t}')
        a('\t\tShading: T')
        a('\t\tCulling: "CullingOff"')
        a('\t}')

    # Materials (one per unique color)
    for hex_color, info in color_to_mat.items():
        r, g, b = hex_to_rgb01(hex_color)
        a(f'\tMaterial: {info["id"]}, "Material::{info["name"]}", "" {{')
        a('\t\tVersion: 102')
        a('\t\tShadingModel: "phong"')
        a('\t\tMultiLayer: 0')
        a('\t\tProperties70:  {')
        a('\t\t\tP: "AmbientColor", "ColorRGB", "Color", "",0,0,0')
        a(f'\t\t\tP: "DiffuseColor", "ColorRGB", "Color", "",{r:.6f},{g:.6f},{b:.6f}')
        a('\t\t\tP: "DiffuseFactor", "Number", "", "A",1')
        a('\t\t\tP: "SpecularColor", "ColorRGB", "Color", "",0.2,0.2,0.2')
        a('\t\t\tP: "SpecularFactor", "Number", "", "A",0.5')
        a('\t\t\tP: "ShininessExponent", "Number", "", "A",20')
        a('\t\t\tP: "Emissive", "Vector3D", "Vector", "",0,0,0')
        a(f'\t\t\tP: "Diffuse", "Vector3D", "Vector", "",{r:.6f},{g:.6f},{b:.6f}')
        a('\t\t}')
        a('\t}')

    a('}')
    a('Connections:  {')
    for p in parts:
        gid = part_ids[id(p)]['gid']
        mid = part_ids[id(p)]['mid']
        mat_id = color_to_mat[p.color_hex]['id']
        name = _safe(p.name)
        a(f'\t;Model::{name}, Model::RootNode')
        a(f'\tC: "OO",{mid},0')
        a(f'\t;Geometry::{name}, Model::{name}')
        a(f'\tC: "OO",{gid},{mid}')
        a(f'\t;Material::{color_to_mat[p.color_hex]["name"]}, Model::{name}')
        a(f'\tC: "OO",{mat_id},{mid}')
    a('}')
    a('Takes:  {')
    a('\tCurrent: ""')
    a('}')
    a('')

    with open(filepath, 'w', encoding='utf-8') as f:
        f.write('\n'.join(out))
    return filepath


# ---------------------------------------------------------------------------
# Public convenience entry point
# ---------------------------------------------------------------------------

def save_3d_as_fbx(mol: Any,
                   filepath: str,
                   atom_color: dict,
                   radius_multi: dict,
                   atom_radius: float = 0.5,
                   pos_multi: float = 1.0,
                   resolution: int = 24,
                   remove_H: bool = True,
                   rotation: Tuple[float, float, float] = (0, 0, 0),
                   bond_color: str | None = None,
                   up_axis: str = 'Z') -> str:
    """Build sphere/cylinder meshes for *mol* and write them as ASCII FBX.

    Returns the final filepath written.
    """
    parts = build_molecule_meshes(
        mol,
        atom_color=atom_color,
        radius_multi=radius_multi,
        atom_radius=atom_radius,
        pos_multi=pos_multi,
        resolution=resolution,
        remove_H=remove_H,
        rotation=rotation,
        bond_color=bond_color,
    )
    return write_ascii_fbx(parts, filepath, up_axis=up_axis)
