"""Build downloadable artefacts (FBX, OBJ zip, interactive HTML, settings).

All exporters take an :class:`Artifacts` request describing the molecule and
the user's appearance choices, build everything in a private temp directory,
and return ``bytes`` ready to feed into ``st.download_button``.
"""
from __future__ import annotations

import io
import logging
import os
import re
import tempfile
import zipfile
from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

from . import config

_LOGGER = logging.getLogger('molecule_finder.exporters')


_SAFE_NAME = re.compile(r'[^A-Za-z0-9._-]+')


def safe_filename(stem: str, fallback: str = 'molecule') -> str:
    """Make *stem* safe for use as a download filename across OSes."""
    cleaned = _SAFE_NAME.sub('_', stem.strip())
    cleaned = cleaned.strip('._-')
    return cleaned or fallback


@dataclass(frozen=True)
class AppearanceSettings:
    palette: Mapping[str, str]                # element symbol → hex color
    resize: Mapping[str, float]               # element symbol → radius multiplier
    atom_radius: float = config.DEFAULT_ATOM_RADIUS
    bond_ratio: float = config.DEFAULT_BOND_RATIO     # bond radius = atom_radius * bond_ratio
    resolution: int = config.DEFAULT_RESOLUTION
    remove_h: bool = True


# ---------------------------------------------------------------------------
# Exporter functions
# ---------------------------------------------------------------------------

def _build_mesh_parts(mol: Any, settings: AppearanceSettings) -> Sequence[Any]:
    from mesh_export import build_molecule_meshes
    return build_molecule_meshes(
        mol,
        atom_color=dict(settings.palette),
        radius_multi=dict(settings.resize),
        atom_radius=settings.atom_radius,
        pos_multi=1.0,
        resolution=max(8, settings.resolution),
        remove_H=settings.remove_h,
    )


def export_fbx(mol: Any, settings: AppearanceSettings, stem: str) -> bytes:
    from mesh_export import write_ascii_fbx
    parts = _build_mesh_parts(mol, settings)
    with tempfile.TemporaryDirectory() as tmp:
        path = write_ascii_fbx(parts, os.path.join(tmp, safe_filename(stem) + '.fbx'))
        with open(path, 'rb') as f:
            data = f.read()
    _LOGGER.info('exported FBX (%d bytes)', len(data))
    return data


def export_obj_zip(mol: Any, settings: AppearanceSettings, stem: str) -> bytes:
    from mesh_export import write_obj
    parts = _build_mesh_parts(mol, settings)
    stem_safe = safe_filename(stem)
    with tempfile.TemporaryDirectory() as tmp:
        obj_path = write_obj(parts, os.path.join(tmp, stem_safe))
        mtl_path = os.path.splitext(obj_path)[0] + '.mtl'
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, 'w', zipfile.ZIP_DEFLATED) as zf:
            zf.write(obj_path, arcname=os.path.basename(obj_path))
            if os.path.exists(mtl_path):
                zf.write(mtl_path, arcname=os.path.basename(mtl_path))
        data = buf.getvalue()
    _LOGGER.info('exported OBJ zip (%d bytes)', len(data))
    return data


def export_html(fig, stem: str) -> bytes:
    """Serialize a plotly figure to a downloadable interactive HTML."""
    cfg = {
        'displaylogo': False,
        'toImageButtonOptions': {
            'format': 'png',
            'filename': safe_filename(stem),
            'scale': 2,
        },
    }
    html = fig.to_html(include_plotlyjs='cdn', full_html=True, config=cfg)
    return html.encode('utf-8')


def export_settings_json(palette: Mapping[str, str],
                         resize: Mapping[str, float],
                         remove_h: bool) -> bytes:
    """Build a settings blob compatible with streamlit-app.py's loader."""
    import json
    blob = {
        'remove_shadow': True,
        'dimension_type': '3D interactive',
        'color_dict': dict(palette),
        'emoji_dict': {},
        'size_multi_slider': 300,
        'resize_dict': dict(resize),
        'update_mol': False,
        'atom_color_select': 'All atoms',
        'use_emoji': False,
        'show_rdkit': False,
        'outline_slider': 1 / 3,
        'reset_size': False,
        'single_bonds': False,
        'change_color_check': True,
        'switch_conf': True,
        'img_format': 'svg',
        'reset_color': False,
        'upload_setting': False,
        'change_size_check': False,
        'removeH': remove_h,
    }
    return json.dumps(blob, indent=2).encode('utf-8')
