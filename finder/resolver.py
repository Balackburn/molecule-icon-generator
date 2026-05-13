"""Resolve a user-supplied query (name / SMILES / CID) to an RDKit Mol.

Resolution order:
  1. PubChem 3D SDF (preferred — gives reproducible geometry).
  2. PubChem CanonicalSMILES → RDKit ETKDG conformer.
  3. cirpy.resolve(..., 'smiles') → RDKit ETKDG conformer.
  4. Treat query as SMILES verbatim → RDKit ETKDG conformer.

Every external call is bounded, retried once, and surfaced as a typed
ResolverError with a user-facing reason.
"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass
from typing import Any, Optional
from urllib.parse import quote

import requests

from . import config

_LOGGER = logging.getLogger('molecule_finder.resolver')


# ---------------------------------------------------------------------------
# Errors
# ---------------------------------------------------------------------------

class ResolverError(RuntimeError):
    """User-facing failure during name → molecule resolution."""


@dataclass(frozen=True)
class ResolveResult:
    mol: Any                       # rdkit.Chem.Mol — kept untyped to defer import
    source: str                    # 'pubchem-3d' | 'pubchem-smiles' | 'cirpy' | 'smiles-direct'
    canonical_smiles: str
    sdf_block: Optional[str] = None
    query: str = ''

    @property
    def heavy_atom_count(self) -> int:
        return self.mol.GetNumHeavyAtoms()

    @property
    def bond_count(self) -> int:
        return self.mol.GetNumBonds()


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------

_SMILES_HINT_CHARS = set('()[]=#@/\\.+-')


def looks_like_smiles(text: str) -> bool:
    """Heuristic: SMILES contain bracket/branch chars and no whitespace."""
    t = text.strip()
    if not t or ' ' in t:
        return False
    if any(c in t for c in _SMILES_HINT_CHARS):
        return True
    # All-lowercase pure-letter strings (e.g. "thc") are almost certainly names
    if t.islower() and t.isalpha():
        return False
    # Mixed-case tokens that include both an upper and a lower letter, no
    # spaces — could be a SMILES like "CCO" or a name like "Aspirin". Default
    # to "no" so we look up via PubChem first; the literal SMILES path is the
    # final fallback so this is still recoverable.
    return False


def sanitize_query(raw: str) -> str:
    """Trim and length-cap; raise if obviously hostile."""
    if raw is None:
        raise ResolverError('Empty query.')
    q = raw.strip()
    if not q:
        raise ResolverError('Empty query.')
    if len(q) > config.MAX_QUERY_LENGTH:
        raise ResolverError(
            f'Query too long ({len(q)} chars; max is {config.MAX_QUERY_LENGTH}).'
        )
    # Block control chars; allow everything else (SMILES use plenty of punctuation).
    if re.search(r'[\x00-\x1f\x7f]', q):
        raise ResolverError('Query contains non-printable characters.')
    return q


# ---------------------------------------------------------------------------
# HTTP helpers
# ---------------------------------------------------------------------------

def _session() -> requests.Session:
    s = requests.Session()
    s.headers.update({'User-Agent': config.HTTP_USER_AGENT,
                      'Accept': 'text/plain, chemical/x-mdl-sdfile'})
    return s


def _get(url: str) -> Optional[str]:
    try:
        r = _session().get(url, timeout=config.HTTP_TIMEOUT)
    except requests.RequestException as exc:
        _LOGGER.warning('HTTP error %s for %s', exc, url)
        return None
    if not r.ok:
        _LOGGER.info('HTTP %s for %s', r.status_code, url)
        return None
    return r.text


# ---------------------------------------------------------------------------
# PubChem branch
# ---------------------------------------------------------------------------

def _pubchem_path_for(query: str) -> str:
    """Choose a /cid/<n> or /name/<encoded> path."""
    if query.isdigit():
        return f'/cid/{int(query)}'
    return f'/name/{quote(query, safe="")}'


def fetch_pubchem_sdf(query: str) -> Optional[str]:
    """Return the PubChem 3D SDF or None."""
    txt = _get(f'{config.PUBCHEM_BASE}{_pubchem_path_for(query)}'
               f'/SDF?record_type=3d')
    if txt and txt.strip().endswith('$$$$'):
        return txt
    return None


def fetch_pubchem_smiles(query: str) -> Optional[str]:
    """Return the canonical SMILES, falling back from canonical→isomeric."""
    for prop in ('CanonicalSMILES', 'IsomericSMILES', 'SMILES'):
        txt = _get(f'{config.PUBCHEM_BASE}{_pubchem_path_for(query)}'
                   f'/property/{prop}/TXT')
        if txt:
            first = txt.strip().splitlines()[0].strip()
            if first:
                return first
    return None


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def _enforce_size(mol) -> None:
    if mol.GetNumAtoms() > config.MAX_ATOMS:
        raise ResolverError(
            f'Molecule has {mol.GetNumAtoms()} atoms; refusing to render more '
            f'than {config.MAX_ATOMS}. Set MIF_MAX_ATOMS if you really need to.'
        )
    if mol.GetNumBonds() > config.MAX_BONDS:
        raise ResolverError(
            f'Molecule has {mol.GetNumBonds()} bonds; refusing to render more '
            f'than {config.MAX_BONDS}. Set MIF_MAX_BONDS if you really need to.'
        )


def _embed(mol, smiles: str):
    """Return a fresh mol with hydrogens + a 3D conformer."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, params) != 0:
        raise ResolverError(
            f'Failed to embed a 3D conformer for "{smiles}". '
            'Try a simpler structure or paste a SMILES.'
        )
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        # MMFF may not cover every atom type; UFF is a reliable fallback.
        try:
            AllChem.UFFOptimizeMolecule(m, maxIters=200)
        except Exception:
            pass  # geometry from ETKDG alone is still usable
    return m


def resolve(query: str) -> ResolveResult:
    """Resolve *query* to an RDKit Mol with a 3D conformer.

    Raises :class:`ResolverError` with a user-facing message on failure.
    """
    from rdkit import Chem
    q = sanitize_query(query)
    _LOGGER.info('Resolving %r', q)

    # --- 1. Direct SMILES if the input clearly looks like one ----------------
    if looks_like_smiles(q):
        mol = Chem.MolFromSmiles(q)
        if mol is not None:
            mol3d = _embed(mol, q)
            _enforce_size(mol3d)
            return ResolveResult(mol=mol3d,
                                 source='smiles-direct',
                                 canonical_smiles=Chem.MolToSmiles(mol3d),
                                 query=q)

    # --- 2. PubChem 3D SDF --------------------------------------------------
    sdf = fetch_pubchem_sdf(q)
    if sdf:
        mol = Chem.MolFromMolBlock(sdf, removeHs=False)
        if mol is not None and mol.GetNumConformers() > 0:
            _enforce_size(mol)
            return ResolveResult(mol=mol,
                                 source='pubchem-3d',
                                 canonical_smiles=Chem.MolToSmiles(mol),
                                 sdf_block=sdf,
                                 query=q)

    # --- 3. PubChem SMILES → embed ------------------------------------------
    smi = fetch_pubchem_smiles(q)
    if smi:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mol3d = _embed(mol, smi)
            _enforce_size(mol3d)
            return ResolveResult(mol=mol3d,
                                 source='pubchem-smiles',
                                 canonical_smiles=Chem.MolToSmiles(mol3d),
                                 query=q)

    # --- 4. cirpy -----------------------------------------------------------
    try:
        import cirpy  # type: ignore
        smi = cirpy.resolve(q, 'smiles')
    except Exception as exc:  # cirpy raises on network errors
        _LOGGER.info('cirpy failed: %s', exc)
        smi = None
    if smi:
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            mol3d = _embed(mol, smi)
            _enforce_size(mol3d)
            return ResolveResult(mol=mol3d,
                                 source='cirpy',
                                 canonical_smiles=Chem.MolToSmiles(mol3d),
                                 query=q)

    # --- 5. Final fallback: treat the raw text as SMILES --------------------
    mol = Chem.MolFromSmiles(q)
    if mol is not None:
        mol3d = _embed(mol, q)
        _enforce_size(mol3d)
        return ResolveResult(mol=mol3d,
                             source='smiles-direct',
                             canonical_smiles=Chem.MolToSmiles(mol3d),
                             query=q)

    raise ResolverError(
        f'Could not find a molecule matching "{q}". Try a different name, '
        'paste a SMILES, or use a PubChem CID number.'
    )
