"""Configuration constants for the Molecule Finder app.

Anything tunable lives here so the UI / library code stays declarative.
Environment variables override the defaults so the same code can run locally,
in Docker, or on a hosted Streamlit deployment.
"""
from __future__ import annotations

import os

# --- HTTP -------------------------------------------------------------------

PUBCHEM_BASE = os.environ.get(
    'MIF_PUBCHEM_BASE',
    'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound',
)
HTTP_TIMEOUT = float(os.environ.get('MIF_HTTP_TIMEOUT', '15'))
HTTP_USER_AGENT = os.environ.get(
    'MIF_USER_AGENT',
    'molecule-icon-generator/1.0 (+https://github.com/Lucandia/molecule-icon-generator)',
)

# --- Validation -------------------------------------------------------------

MAX_QUERY_LENGTH = int(os.environ.get('MIF_MAX_QUERY_LENGTH', '200'))
MAX_ATOMS = int(os.environ.get('MIF_MAX_ATOMS', '600'))
MAX_BONDS = int(os.environ.get('MIF_MAX_BONDS', '600'))

# --- Caching ----------------------------------------------------------------

# Streamlit cache TTL (seconds). PubChem responses are stable, so a day is fine.
CACHE_TTL_SECONDS = int(os.environ.get('MIF_CACHE_TTL', '86400'))

# --- Render defaults --------------------------------------------------------

DEFAULT_RESOLUTION = int(os.environ.get('MIF_DEFAULT_RESOLUTION', '24'))
PREVIEW_RESOLUTION = int(os.environ.get('MIF_PREVIEW_RESOLUTION', '20'))
DEFAULT_ATOM_RADIUS = float(os.environ.get('MIF_DEFAULT_ATOM_RADIUS', '0.5'))
DEFAULT_BOND_RATIO = float(os.environ.get('MIF_DEFAULT_BOND_RATIO', '0.25'))

# --- Presets ----------------------------------------------------------------

PRESETS: dict[str, str] = {
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

# Order matters here — these show in the sidebar as fields.
PALETTE_KEYS: tuple[str, ...] = (
    'C', 'H', 'N', 'O', 'S', 'P',
    'F', 'Cl', 'Br', 'I',
    'Fe', 'Zn', 'Cu', 'Mg', 'Ca', 'Na', 'K',
    'other', 'Bond', 'Background',
)

# --- Logging ----------------------------------------------------------------

LOG_LEVEL = os.environ.get('MIF_LOG_LEVEL', 'INFO').upper()
