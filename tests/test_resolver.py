"""Tests for finder.resolver — hits PubChem live, marked as integration."""
from __future__ import annotations

import os
import sys

import pytest

HERE = os.path.dirname(__file__)
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)

from finder import resolver as r  # noqa: E402


# ---------------------------------------------------------------------------
# Pure-Python helpers
# ---------------------------------------------------------------------------

def test_sanitize_query_strips_whitespace():
    assert r.sanitize_query('  caffeine \n') == 'caffeine'


def test_sanitize_query_rejects_empty():
    with pytest.raises(r.ResolverError):
        r.sanitize_query('')


def test_sanitize_query_rejects_overlong():
    with pytest.raises(r.ResolverError):
        r.sanitize_query('a' * 10_000)


def test_sanitize_query_rejects_control_chars():
    with pytest.raises(r.ResolverError):
        r.sanitize_query('caffeine\x00')


@pytest.mark.parametrize('text,expected', [
    ('CC(=O)OC1=CC=CC=C1C(=O)O', True),    # aspirin SMILES
    ('CCO', False),                          # ambiguous; treated as a name first
    ('caffeine', False),
    ('THC', False),
    ('C[C@H](N)C(=O)O', True),               # alanine
])
def test_looks_like_smiles(text, expected):
    assert r.looks_like_smiles(text) is expected


def test_pubchem_path_for_uses_cid_when_numeric():
    assert r._pubchem_path_for('16078') == '/cid/16078'
    assert r._pubchem_path_for('caffeine').startswith('/name/')


# ---------------------------------------------------------------------------
# Live PubChem integration (skipped if offline / MIF_OFFLINE=1)
# ---------------------------------------------------------------------------

OFFLINE = os.environ.get('MIF_OFFLINE') == '1'
integration = pytest.mark.skipif(OFFLINE,
                                 reason='MIF_OFFLINE=1 — skipping live HTTP test')


@integration
@pytest.mark.parametrize('query', ['caffeine', 'aspirin', '2244'])
def test_resolve_smoke(query):
    pytest.importorskip('rdkit')
    res = r.resolve(query)
    assert res.mol.GetNumConformers() == 1
    assert res.heavy_atom_count > 0
    assert res.canonical_smiles
    assert res.source in {'pubchem-3d', 'pubchem-smiles', 'cirpy',
                          'smiles-direct'}


@integration
def test_resolve_smiles_direct():
    pytest.importorskip('rdkit')
    res = r.resolve('CC(=O)OC1=CC=CC=C1C(=O)O')
    assert res.source == 'smiles-direct'
    assert res.heavy_atom_count == 13


@integration
def test_resolve_unknown_raises():
    pytest.importorskip('rdkit')
    with pytest.raises(r.ResolverError):
        r.resolve('this is not a real molecule name xyzzy fnord')
