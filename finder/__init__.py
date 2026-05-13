"""Production-grade Molecule Finder package.

Contains:
  - finder.config     — runtime configuration + presets
  - finder.resolver   — name/SMILES/CID → RDKit mol (PubChem + cirpy)
  - finder.exporters  — file builders for FBX, OBJ and interactive HTML
  - finder.app        — the Streamlit UI

Run with:
    python -m streamlit run streamlit_finder.py
"""
from __future__ import annotations

__all__ = ['config', 'resolver', 'exporters']
__version__ = '1.0.0'
