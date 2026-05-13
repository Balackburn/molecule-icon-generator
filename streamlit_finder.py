"""Molecule Finder — production-grade Streamlit front-end.

Type a molecule name, SMILES, or PubChem CID. Preview the 3D structure inline
and download as FBX / OBJ+MTL (.zip) / interactive HTML, plus a settings
JSON re-usable with the original streamlit-app.py.

Run with:
    python -m streamlit run streamlit_finder.py

Configuration is environment-driven; see finder/config.py for the keys.
"""
from __future__ import annotations

import json
import logging
import time

import streamlit as st

import molecule_icon_generator as mig
from finder import config
from finder.exporters import (AppearanceSettings, export_fbx, export_html,
                              export_obj_zip, export_settings_json,
                              safe_filename)
from finder.resolver import ResolverError, resolve


# ---------------------------------------------------------------------------
# Logging — once, at import time
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=getattr(logging, config.LOG_LEVEL, logging.INFO),
    format='%(asctime)s %(name)s %(levelname)s %(message)s',
)
_LOGGER = logging.getLogger('molecule_finder.app')


# ---------------------------------------------------------------------------
# Cache-key helpers (JSON is hashable and safe to round-trip)
# ---------------------------------------------------------------------------

def _pack(mapping: dict) -> str:
    return json.dumps(mapping, sort_keys=True)


def _unpack(packed: str) -> dict:
    return json.loads(packed)


# ---------------------------------------------------------------------------
# Cached helpers
# ---------------------------------------------------------------------------

@st.cache_data(ttl=config.CACHE_TTL_SECONDS, show_spinner=False)
def cached_resolve(query: str):
    """Cache the resolve() result by query string."""
    return resolve(query)


@st.cache_data(ttl=config.CACHE_TTL_SECONDS, show_spinner=False)
def cached_figure(canonical_smiles: str, palette_key: str, resize_key: str,
                  resolution: int, remove_h: bool, bg: str):
    """Build the plotly figure. Keyed on canonical SMILES + cosmetic settings."""
    palette = _unpack(palette_key)
    resize = {k: float(v) for k, v in _unpack(resize_key).items()}
    mol = _mol_from_smiles(canonical_smiles)
    fig = mig.graph_3d(
        mol,
        atom_color=palette,
        radius_multi=resize,
        atom_radius=100,
        pos_multi=300,
        resolution=max(8, resolution),
        remove_H=remove_h,
    )
    fig.update_layout(
        height=620,
        paper_bgcolor=bg,
        scene=dict(bgcolor=bg, aspectmode='data'),
        margin=dict(l=0, r=0, t=0, b=0),
    )
    return fig


def _mol_from_smiles(smiles: str):
    """Build a 3D-embedded mol from a SMILES string (with hydrogens)."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise RuntimeError(f'Invalid SMILES: {smiles!r}')
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, params) != 0:
        raise RuntimeError(f'Could not embed a 3D conformer for {smiles!r}.')
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        try:
            AllChem.UFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass
    return m


# ---------------------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------------------

def render_sidebar() -> AppearanceSettings:
    st.sidebar.header('Customization')

    remove_h = st.sidebar.checkbox('Remove hydrogens', value=True,
                                   help='Hide non-chiral H atoms.')
    resolution = st.sidebar.slider('Mesh resolution', 8, 64,
                                   config.DEFAULT_RESOLUTION,
                                   help='Triangulation segments per sphere/cylinder.')
    atom_radius = st.sidebar.slider('Atom radius', 0.1, 1.5,
                                    config.DEFAULT_ATOM_RADIUS, 0.05)
    bond_ratio = st.sidebar.slider('Bond thickness ÷ atom radius',
                                   0.05, 0.5, config.DEFAULT_BOND_RATIO, 0.01)

    st.sidebar.divider()
    st.sidebar.subheader('Palette')

    base = dict(mig.color_map)
    base.setdefault('Bond', '#979797')
    base.setdefault('Background', '#FFFFFF')
    base.setdefault('other', '#FF66CC')

    mode = st.sidebar.radio(
        'Color scheme',
        ['CPK (default)', 'Monochrome', 'Custom'],
        index=0,
    )

    palette: dict[str, str] = dict(base)

    if mode == 'Monochrome':
        mono = st.sidebar.color_picker('All atoms', '#FC65B6')
        palette = {k: mono for k in base.keys()}
        palette['Bond'] = st.sidebar.color_picker('Bond', '#979797')
        palette['Background'] = st.sidebar.color_picker('Background', '#FFFFFF')
    elif mode == 'Custom':
        cols = st.sidebar.columns(2)
        for i, k in enumerate(config.PALETTE_KEYS):
            with cols[i % 2]:
                palette[k] = st.color_picker(k, base.get(k, '#888888'),
                                             key=f'col_{k}')
    else:
        palette['Background'] = '#FFFFFF'

    return AppearanceSettings(
        palette=palette,
        resize=dict(mig.atom_resize),
        atom_radius=atom_radius,
        bond_ratio=bond_ratio,
        resolution=resolution,
        remove_h=remove_h,
    )


# ---------------------------------------------------------------------------
# Header + search
# ---------------------------------------------------------------------------

def render_header() -> str:
    st.title('🧪 Molecule Finder · 3D preview & download')
    st.caption(
        'Type a molecule name (e.g. **THC**, **caffeine**), paste a SMILES, '
        'or enter a PubChem CID. The 3D structure renders inline and you can '
        'download it as FBX, OBJ + MTL, interactive HTML, or settings JSON.'
    )

    if 'query' not in st.session_state:
        st.session_state['query'] = 'tetrahydrocannabinol'

    with st.form('search_form', clear_on_submit=False, border=False):
        col_q, col_btn = st.columns([6, 1])
        with col_q:
            new_q = st.text_input(
                'Molecule name, SMILES or PubChem CID',
                value=st.session_state['query'],
                max_chars=config.MAX_QUERY_LENGTH,
                label_visibility='collapsed',
                placeholder='Try "THC", "caffeine" or "CC(=O)OC1=CC=CC=C1C(=O)O"',
            )
        with col_btn:
            submitted = st.form_submit_button('Find', use_container_width=True,
                                              type='primary')
    if submitted and new_q.strip():
        st.session_state['query'] = new_q.strip()

    st.caption('Quick picks:')
    preset_cols = st.columns(len(config.PRESETS))
    for i, (label, term) in enumerate(config.PRESETS.items()):
        with preset_cols[i]:
            if st.button(label, key=f'preset_{label}', use_container_width=True):
                st.session_state['query'] = term
                st.rerun()

    return st.session_state['query']


# ---------------------------------------------------------------------------
# Main rendering
# ---------------------------------------------------------------------------

def render_result(query: str, settings: AppearanceSettings) -> None:
    t0 = time.perf_counter()
    try:
        with st.spinner(f'Looking up "{query}"…'):
            result = cached_resolve(query)
    except ResolverError as exc:
        st.error(str(exc), icon='⚠️')
        st.stop()
    except Exception as exc:  # last-resort safety net
        _LOGGER.exception('unexpected resolver failure')
        st.error(
            f'Unexpected resolver error: {exc}. '
            'Please retry or report this with the query you used.',
            icon='⚠️',
        )
        st.stop()

    resolve_secs = time.perf_counter() - t0

    # --- Header strip with quick facts -------------------------------------
    c1, c2, c3, c4 = st.columns(4)
    c1.metric('Heavy atoms', result.heavy_atom_count)
    c2.metric('Bonds', result.bond_count)
    c3.metric('Source', result.source)
    c4.metric('Resolved in', f'{resolve_secs:0.2f} s')

    with st.expander('Resolved details', expanded=False):
        st.markdown('**Canonical SMILES**')
        st.code(result.canonical_smiles, language='text')
        if result.sdf_block:
            st.markdown('**PubChem 3D SDF**')
            st.text_area('SDF', result.sdf_block, height=180,
                         label_visibility='collapsed')

    # --- 3D preview --------------------------------------------------------
    with st.spinner('Rendering 3D mesh…'):
        try:
            fig = cached_figure(
                canonical_smiles=result.canonical_smiles,
                palette_key=_pack(dict(settings.palette)),
                resize_key=_pack(dict(settings.resize)),
                resolution=config.PREVIEW_RESOLUTION,
                remove_h=settings.remove_h,
                bg=settings.palette.get('Background', '#FFFFFF'),
            )
        except Exception as exc:
            _LOGGER.exception('plotly figure build failed')
            st.error(f'Failed to render preview: {exc}', icon='⚠️')
            st.stop()

    preview_col, dl_col = st.columns([3, 1])

    with preview_col:
        st.plotly_chart(
            fig,
            use_container_width=True,
            config={'displaylogo': False,
                    'toImageButtonOptions': {
                        'format': 'png',
                        'filename': safe_filename(query),
                        'scale': 2,
                    }},
        )

    with dl_col:
        st.subheader('Downloads')
        stem = safe_filename(query)

        def _safe_download(label, builder, file_name, mime):
            try:
                data = builder()
            except Exception as exc:
                _LOGGER.exception('export failed: %s', label)
                st.warning(f'{label}: {exc}', icon='⚠️')
                return
            st.download_button(label, data=data, file_name=file_name,
                               mime=mime, use_container_width=True)

        _safe_download(
            '⬇︎ FBX (3D model)',
            lambda: export_fbx(result.mol, settings, stem),
            file_name=f'{stem}.fbx',
            mime='application/octet-stream',
        )
        _safe_download(
            '⬇︎ OBJ + MTL (.zip)',
            lambda: export_obj_zip(result.mol, settings, stem),
            file_name=f'{stem}.zip',
            mime='application/zip',
        )
        _safe_download(
            '⬇︎ Interactive HTML',
            lambda: export_html(fig, stem),
            file_name=f'{stem}.html',
            mime='text/html',
        )
        _safe_download(
            '⬇︎ Settings JSON',
            lambda: export_settings_json(settings.palette, settings.resize,
                                         settings.remove_h),
            file_name='molecule_icon_settings.json',
            mime='application/json',
        )


def render_footer():
    st.divider()
    st.caption(
        'Built on top of [Lucandia/Molecule-Icon-Generator]'
        '(https://github.com/Lucandia/molecule-icon-generator). '
        'Names resolve via PubChem REST with a cirpy + direct-SMILES fallback. '
        'FBX export is pure Python — no Autodesk SDK required.'
    )


def main():
    st.set_page_config(
        page_title='Molecule Finder · 3D preview & FBX export',
        page_icon='🧪',
        layout='wide',
        menu_items={
            'Get help': 'https://github.com/Lucandia/molecule-icon-generator',
            'Report a bug': 'https://github.com/Lucandia/molecule-icon-generator/issues',
            'About': 'Molecule Finder — type a name, get a 3D model.',
        },
    )
    query = render_header()
    settings = render_sidebar()
    render_result(query, settings)
    render_footer()


main()
