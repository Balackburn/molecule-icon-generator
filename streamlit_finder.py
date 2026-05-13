"""Molecule Finder — type a molecule name, preview it in 3D, download.

A streamlined companion to streamlit-app.py focused on:
  1. Name → SMILES resolution (cirpy → PubChem fallback)
  2. Inline 3D preview (plotly)
  3. One-click downloads (FBX, OBJ, interactive HTML, snapshot PNG/SVG)
  4. Lightweight customization (atom palette + bond color + size + resolution)

Run with:  python -m streamlit run streamlit_finder.py
"""

from __future__ import annotations

import json
import os
import tempfile
from io import BytesIO

import requests
import streamlit as st

import molecule_icon_generator as mig
import mesh_export as mx


# ---------------------------------------------------------------------------
# Constants / presets
# ---------------------------------------------------------------------------

PUBCHEM_REST = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound'

PRESETS = {
    'THC':         'tetrahydrocannabinol',
    'CBD':         'cannabidiol',
    'Caffeine':    'caffeine',
    'Aspirin':     'aspirin',
    'Paracetamol': 'paracetamol',
    'Glucose':     'glucose',
    'Adrenaline':  'adrenaline',
    'Serotonin':   'serotonin',
    'Dopamine':    'dopamine',
    'Nicotine':    'nicotine',
    'Cholesterol': 'cholesterol',
    'Vitamin C':   'ascorbic acid',
}

# Common elements + an "other" bucket for everything else
PALETTE_KEYS = ['C', 'H', 'N', 'O', 'S', 'P',
                'F', 'Cl', 'Br', 'I',
                'Fe', 'Zn', 'Cu', 'Mg', 'Ca', 'Na', 'K',
                'other', 'Bond', 'Background']


# ---------------------------------------------------------------------------
# Resolvers — name to SMILES / 3D SDF
# ---------------------------------------------------------------------------

@st.cache_data(show_spinner=False)
def resolve_name_to_smiles(name: str) -> str | None:
    """Try cirpy first, then PubChem REST. Returns SMILES string or None."""
    name = name.strip()
    if not name:
        return None
    # SMILES heuristic: if string contains typical SMILES characters and no spaces,
    # accept it as-is.
    if any(c in name for c in '()[]=#') and ' ' not in name:
        return name
    try:
        import cirpy  # type: ignore
        s = cirpy.resolve(name, 'smiles')
        if s:
            return s
    except Exception:
        pass
    # PubChem fallback
    try:
        url = f'{PUBCHEM_REST}/name/{requests.utils.quote(name)}/property/CanonicalSMILES/TXT'
        r = requests.get(url, timeout=15)
        if r.ok:
            return r.text.strip().splitlines()[0]
    except Exception:
        pass
    return None


@st.cache_data(show_spinner=False)
def fetch_pubchem_sdf(name: str) -> str | None:
    """Fetch a 3D SDF from PubChem for a name. Returns SDF text or None."""
    try:
        url = f'{PUBCHEM_REST}/name/{requests.utils.quote(name)}/SDF?record_type=3d'
        r = requests.get(url, timeout=20)
        if r.ok and r.text.strip().endswith('$$$$'):
            return r.text
    except Exception:
        pass
    return None


def mol_from_query(query: str):
    """Resolve query (name or SMILES) → RDKit mol with a 3D conformer."""
    sdf = fetch_pubchem_sdf(query)
    if sdf:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        m = Chem.MolFromMolBlock(sdf)
        if m is not None:
            # PubChem 3D SDFs already include H + 3D coords — keep them.
            return m, query, sdf
    smiles = resolve_name_to_smiles(query)
    if not smiles:
        return None, None, None
    mol = mig.parse_structure(smiles, nice_conformation=False, dimension_3=True)
    return mol, smiles, None


# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------

st.set_page_config(page_title='Molecule Finder (3D + FBX)', layout='wide')

st.title('🧪 Molecule Finder · 3D preview & download')
st.caption(
    'Type a molecule name (e.g. THC, caffeine) or paste a SMILES. '
    'Preview the 3D structure and download as FBX / OBJ / HTML / snapshot.'
)


# ---------------------------------------------------------------------------
# Search row + presets
# ---------------------------------------------------------------------------

if 'query' not in st.session_state:
    st.session_state['query'] = 'tetrahydrocannabinol'

col_search, col_btn = st.columns([6, 1])
with col_search:
    query = st.text_input('Molecule name or SMILES',
                          value=st.session_state['query'],
                          key='query_input',
                          help='Try "THC", "caffeine", "aspirin" or a SMILES string.')
with col_btn:
    st.write('')  # spacer
    st.write('')
    go = st.button('Find', use_container_width=True, type='primary')

st.caption('Quick picks:')
preset_cols = st.columns(len(PRESETS))
for i, (label, term) in enumerate(PRESETS.items()):
    with preset_cols[i]:
        if st.button(label, key=f'preset_{label}', use_container_width=True):
            st.session_state['query'] = term
            st.rerun()

if go:
    st.session_state['query'] = query

active_query = st.session_state['query']


# ---------------------------------------------------------------------------
# Sidebar — customization
# ---------------------------------------------------------------------------

st.sidebar.header('Customization')

remove_H = st.sidebar.checkbox('Remove hydrogens', value=True)
resolution = st.sidebar.slider('Mesh resolution', 8, 64, 28,
                               help='Triangulation segments for spheres and cylinders.')
atom_radius = st.sidebar.slider('Atom radius', 0.1, 1.5, 0.5, 0.05)
bond_radius_factor = st.sidebar.slider('Bond thickness × atom radius',
                                       0.05, 0.5, 0.25, 0.01)

st.sidebar.divider()
st.sidebar.subheader('Colors')

default_colors = mig.color_map
defaults = {k: default_colors.get(k, '#888888') for k in PALETTE_KEYS}
defaults['Bond'] = default_colors.get('Bond', '#979797')
defaults['Background'] = '#FFFFFF'
defaults['other'] = '#FF66CC'

palette_mode = st.sidebar.radio(
    'Palette',
    ['CPK (default)', 'Monochrome', 'Custom'],
    horizontal=False,
)

palette: dict[str, str] = dict(default_colors)

if palette_mode == 'Monochrome':
    mono = st.sidebar.color_picker('All atoms', '#FC65B6')
    palette = {k: mono for k in default_colors.keys()}
    palette['Bond'] = st.sidebar.color_picker('Bond', '#979797')
    palette['Background'] = st.sidebar.color_picker('Background', '#FFFFFF')
elif palette_mode == 'Custom':
    cols_left = st.sidebar.columns(2)
    palette = dict(default_colors)
    for i, k in enumerate(PALETTE_KEYS):
        with cols_left[i % 2]:
            palette[k] = st.color_picker(k, defaults.get(k, '#888888'),
                                         key=f'col_{k}')
else:
    palette['Background'] = '#FFFFFF'


# ---------------------------------------------------------------------------
# Resolve + build mol
# ---------------------------------------------------------------------------

mol, resolved, sdf_blob = (None, None, None)
err = None
with st.spinner(f'Looking up "{active_query}"…'):
    try:
        mol, resolved, sdf_blob = mol_from_query(active_query)
    except Exception as e:
        err = str(e)

if err:
    st.error(f'Failed to resolve "{active_query}": {err}')
    st.stop()
if mol is None:
    st.error(
        f'Could not find a molecule for "{active_query}".  Try a different '
        'name or paste a SMILES string.'
    )
    st.stop()

# Cards: resolved info
info_a, info_b, info_c = st.columns(3)
info_a.metric('Atoms (after H rule)',
              mol.GetNumHeavyAtoms() if remove_H else mol.GetNumAtoms())
info_b.metric('Bonds', mol.GetNumBonds())
info_c.metric('Source', 'PubChem 3D' if sdf_blob else 'RDKit ETKDG')

with st.expander('Resolved SMILES / Mol block', expanded=False):
    from rdkit import Chem
    st.code(Chem.MolToSmiles(mol), language='text')
    if sdf_blob:
        st.text_area('PubChem SDF', sdf_blob, height=180)


# ---------------------------------------------------------------------------
# Build meshes (used for both download and inline preview)
# ---------------------------------------------------------------------------

resize_dict = dict(mig.atom_resize)

with st.spinner('Building 3D meshes…'):
    fig = mig.graph_3d(
        mol,
        atom_color=palette,
        radius_multi=resize_dict,
        atom_radius=100,            # mig graph_3d uses pixel-ish units
        pos_multi=300,
        resolution=max(8, resolution),
        remove_H=remove_H,
    )
    # tighten the preview camera + background
    fig.update_layout(
        height=620,
        paper_bgcolor=palette['Background'],
        scene=dict(bgcolor=palette['Background']),
        margin=dict(l=0, r=0, t=0, b=0),
    )


# ---------------------------------------------------------------------------
# Preview + downloads
# ---------------------------------------------------------------------------

preview_col, download_col = st.columns([3, 1])

with preview_col:
    st.plotly_chart(fig, use_container_width=True,
                    config={'displaylogo': False,
                            'toImageButtonOptions': {
                                'format': 'png',
                                'filename': active_query.replace(' ', '_'),
                                'scale': 2,
                            }})

with download_col:
    st.subheader('Downloads')

    base = os.path.join(tempfile.gettempdir(),
                        active_query.replace(' ', '_') or 'molecule')

    # FBX
    try:
        fbx_path = mig.save_3d_fbx(
            mol,
            name=os.path.basename(base),
            directory=os.path.dirname(base) or '.',
            atom_color=palette,
            radius_multi=resize_dict,
            atom_radius=atom_radius,
            pos_multi=1.0,
            resolution=max(8, resolution),
            remove_H=remove_H,
        )
        with open(fbx_path, 'rb') as f:
            st.download_button('⬇︎ FBX (3D model)',
                               data=f,
                               file_name=os.path.basename(fbx_path),
                               mime='application/octet-stream',
                               use_container_width=True)
    except Exception as e:
        st.warning(f'FBX failed: {e}')

    # OBJ + MTL bundle
    try:
        from mesh_export import build_molecule_meshes, write_obj
        parts = build_molecule_meshes(
            mol,
            atom_color=palette,
            radius_multi=resize_dict,
            atom_radius=atom_radius,
            pos_multi=1.0,
            resolution=max(8, resolution),
            remove_H=remove_H,
        )
        obj_base = base
        obj_path = write_obj(parts, obj_base)
        mtl_path = os.path.splitext(obj_path)[0] + '.mtl'

        # Bundle .obj + .mtl into a single .zip for the download button
        import zipfile
        zip_buf = BytesIO()
        with zipfile.ZipFile(zip_buf, 'w', zipfile.ZIP_DEFLATED) as zf:
            zf.write(obj_path, arcname=os.path.basename(obj_path))
            zf.write(mtl_path, arcname=os.path.basename(mtl_path))
        st.download_button('⬇︎ OBJ + MTL (.zip)',
                           data=zip_buf.getvalue(),
                           file_name=os.path.basename(base) + '.zip',
                           mime='application/zip',
                           use_container_width=True)
    except Exception as e:
        st.warning(f'OBJ failed: {e}')

    # Interactive HTML
    try:
        html_path = base + '.html'
        fig.write_html(html_path,
                       config={'toImageButtonOptions':
                               {'format': 'png',
                                'filename': os.path.basename(base),
                                'scale': 2}})
        with open(html_path, 'rb') as f:
            st.download_button('⬇︎ Interactive HTML',
                               data=f,
                               file_name=os.path.basename(html_path),
                               mime='text/html',
                               use_container_width=True)
    except Exception as e:
        st.warning(f'HTML failed: {e}')

    # Settings JSON — re-usable with the original streamlit-app
    settings_json = {
        'remove_shadow': True,
        'dimension_type': '3D interactive',
        'color_dict': palette,
        'emoji_dict': {},
        'size_multi_slider': 300,
        'resize_dict': resize_dict,
        'update_mol': False,
        'atom_color_select': 'All atoms',
        'use_emoji': False,
        'show_rdkit': False,
        'outline_slider': 1/3,
        'reset_size': False,
        'single_bonds': False,
        'change_color_check': True,
        'switch_conf': True,
        'img_format': 'svg',
        'reset_color': False,
        'upload_setting': False,
        'change_size_check': False,
        'removeH': remove_H,
    }
    st.download_button('⬇︎ Settings JSON',
                       data=json.dumps(settings_json, indent=2),
                       file_name='molecule_icon_settings.json',
                       mime='application/json',
                       use_container_width=True)


st.divider()
st.caption(
    'Built on top of [Lucandia/Molecule-Icon-Generator]'
    '(https://github.com/Lucandia/molecule-icon-generator). '
    'FBX export written in pure Python; resolves names via cirpy '
    'and the PubChem REST API.'
)
