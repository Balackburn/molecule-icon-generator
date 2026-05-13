"""Unit tests for mesh_export — pure-Python sphere/cylinder/FBX/OBJ writers.

These tests do not require RDKit.
"""
from __future__ import annotations

import os
import re
import sys
import tempfile

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

import mesh_export as mx  # noqa: E402


# ---------------------------------------------------------------------------
# Primitive builders
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('segments,rings_expected_lo', [(8, 4), (16, 8), (32, 16)])
def test_sphere_mesh_topology(segments, rings_expected_lo):
    v, n, t = mx.sphere_mesh((1, 2, 3), radius=0.5, segments=segments)
    # All normals are unit length
    norms = np.linalg.norm(n, axis=1)
    assert np.allclose(norms, 1.0, atol=1e-5)
    # Every triangle index is in range
    assert t.min() >= 0
    assert t.max() < v.shape[0]
    # Centered roughly at (1,2,3)
    assert np.allclose(v.mean(axis=0), (1.0, 2.0, 3.0), atol=0.1)


def test_cylinder_zero_length_is_empty():
    v, n, t = mx.cylinder_mesh((0, 0, 0), (0, 0, 0), radius=0.1)
    assert v.shape == (0, 3)
    assert t.shape == (0, 3)


def test_cylinder_mesh_has_caps():
    v, n, t = mx.cylinder_mesh((0, 0, 0), (1, 0, 0), radius=0.2, segments=16)
    # Side ring vertices have normals perpendicular to axis (x-axis); cap
    # vertices have normals along ±x.
    along_x = np.abs(n[:, 0])
    # at least some normals point along x (the caps)
    assert (along_x > 0.99).sum() > 16


# ---------------------------------------------------------------------------
# FBX writer
# ---------------------------------------------------------------------------

def _three_parts():
    a_v, a_n, a_t = mx.sphere_mesh((0, 0, 0), 0.5, segments=12)
    b_v, b_n, b_t = mx.sphere_mesh((1, 0, 0), 0.4, segments=12)
    c_v, c_n, c_t = mx.cylinder_mesh((0, 0, 0), (1, 0, 0), radius=0.05,
                                     segments=12)
    return [
        mx.MeshPart('atom_C0', a_v, a_n, a_t, '#222222'),
        mx.MeshPart('atom_O1', b_v, b_n, b_t, '#FF0000'),
        mx.MeshPart('bond_0',  c_v, c_n, c_t, '#888888'),
    ]


def test_fbx_writer_produces_balanced_consistent_file(tmp_path):
    parts = _three_parts()
    out = tmp_path / 'mol.fbx'
    mx.write_ascii_fbx(parts, str(out))
    text = out.read_text()

    assert text.startswith('; FBX 7.4.0')
    assert text.count('{') == text.count('}')

    for label in ('Vertices', 'PolygonVertexIndex', 'Normals'):
        for m in re.finditer(rf'{label}: \*(\d+) \{{\s*a:\s*([^}}]+)\}}', text):
            declared = int(m.group(1))
            actual = len(m.group(2).strip().split(','))
            assert declared == actual, f'{label}: {declared} != {actual}'

    # 3 unique colors → 3 materials
    assert len(re.findall(r'^\tMaterial: ', text, re.MULTILINE)) == 3
    # 3 parts × 3 OO connections each
    assert len(re.findall(r'^\tC: "OO",', text, re.MULTILINE)) == 9


def test_fbx_writer_dedups_materials_by_color(tmp_path):
    parts = _three_parts()
    # Force first two atoms to share a color
    parts[0] = mx.MeshPart('atom_a', parts[0].vertices, parts[0].normals,
                           parts[0].triangles, '#FF0000')
    parts[1] = mx.MeshPart('atom_b', parts[1].vertices, parts[1].normals,
                           parts[1].triangles, '#FF0000')
    out = tmp_path / 'shared.fbx'
    mx.write_ascii_fbx(parts, str(out))
    text = out.read_text()
    assert len(re.findall(r'^\tMaterial: ', text, re.MULTILINE)) == 2


def test_fbx_appends_extension_when_missing(tmp_path):
    parts = _three_parts()
    out = tmp_path / 'no_extension'
    written = mx.write_ascii_fbx(parts, str(out))
    assert written.endswith('.fbx')


def test_fbx_writer_rejects_empty_part_list(tmp_path):
    with pytest.raises(ValueError):
        mx.write_ascii_fbx([], str(tmp_path / 'x.fbx'))


# ---------------------------------------------------------------------------
# OBJ writer
# ---------------------------------------------------------------------------

def test_obj_writer_counts(tmp_path):
    parts = _three_parts()
    obj_path = mx.write_obj(parts, str(tmp_path / 'mol'))
    text = open(obj_path).read()
    expected_v = sum(p.vertices.shape[0] for p in parts)
    expected_f = sum(p.triangles.shape[0] for p in parts)
    assert text.count('\nv ') == expected_v
    assert text.count('\nf ') == expected_f
    # MTL exists and references both materials
    mtl_path = os.path.splitext(obj_path)[0] + '.mtl'
    mtl = open(mtl_path).read()
    assert mtl.count('newmtl ') == 3  # one per distinct color


# ---------------------------------------------------------------------------
# Color helper
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('hex_in,expected', [
    ('#000000', (0.0, 0.0, 0.0)),
    ('#FFFFFF', (1.0, 1.0, 1.0)),
    ('#fc65b6', (252/255, 101/255, 182/255)),
    ('bad-color', (0.5, 0.5, 0.5)),
    (None, (0.5, 0.5, 0.5)),
])
def test_hex_to_rgb01(hex_in, expected):
    got = mx.hex_to_rgb01(hex_in)
    for g, e in zip(got, expected):
        assert abs(g - e) < 1e-6
