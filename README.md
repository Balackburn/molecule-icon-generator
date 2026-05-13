# Molecule icon Generator!
<img src="example/paracetamol.jpeg" width=300 height=300>
Generate nice icons of molecules from SMILES.

This program follows the topology of SMILES chemical structures and creates an icon.
The atoms' colours are inspired by the [CPK colouring convention](https://sciencenotes.org/molecule-atom-colors-cpk-colors/)

## Try the web app:

[Molecules-icons generator web app](https://molecule-icon-generator.streamlit.app/) powered by streamlit

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://molecule-icon-generator.streamlit.app/)

## Cite this work:

If you use the icons please cite this project:

```
Lucandia. Lucandia/Molecule-Icon-Generator; Zenodo, 2022. https://doi.org/10.5281/ZENODO.7388429.
```
[![DOI](https://zenodo.org/badge/530035520.svg)](https://zenodo.org/badge/latestdoi/530035520)

## Run it on local:

### Step 1: clone the repository

```
git clone https://github.com/lmonari5/molecule-icon-generator.git
```

### Step 2: install packages and requirements

For Linux:

```
xargs -a packages.txt sudo apt-get install 
pip install -r requirements.txt
```

### Step 3: Have Fun

- Run the code from the command line:

 ```
  python molecule_icon_generator.py "CC(=O)Nc1ccc(cc1)O" --name paracetamol --rdkit_svg
 ```

- Or from the python interpreter:

 ```
 import molecule_icon_generator as mig 
 molecule = mig.parse_structure("CC(=O)Nc1ccc(cc1)O", remove_H=False)
 mig.icon_print(molecule, name = 'paracetamol', rdkit_svg = True, single_bonds = False, remove_H = False, verbose=False)
 ```

- Or use the Streamlit functionalities! From the terminal run:

 ```
python -m streamlit run streamlit_app.py
 ```

## Export to 3D (FBX / OBJ)

Save the molecule as a real 3D model with sphere atoms and cylinder bonds.
The export lives in `mesh_export.py` and only uses `numpy`, `scipy` and
`rdkit` — already in `requirements.txt`. ASCII FBX 7.4 is written so the file
opens in Blender, Maya, 3ds Max, Cinema 4D, Unity and Unreal without needing
the Autodesk FBX SDK. Materials are deduplicated by atom color, so all
carbons share one material, all hydrogens share another, etc.

Command line:

```
# write paracetamol.fbx (with a 3D conformer) next to the 2D icon
python molecule_icon_generator.py "CC(=O)Nc1ccc(cc1)O" --name paracetamol --fbx

# also write a Wavefront OBJ + MTL alongside the FBX, and bump tessellation
python molecule_icon_generator.py "CCO" --name ethanol --fbx --obj --fbx_resolution 32
```

From Python:

```
import molecule_icon_generator as mig
mol = mig.parse_structure("CC(=O)Nc1ccc(cc1)O", dimension_3=True)
mig.save_3d_fbx(mol, name="paracetamol", directory=".", verbose=True)
```

In the Streamlit app, switch the dimension selector to `3D interactive` and
a `Download 3D model (FBX)` button appears next to the existing HTML
download.

## Molecule Finder (production-grade companion app)

A streamlined Streamlit front-end: type a molecule name (e.g. `THC`,
`caffeine`), a SMILES string, or a PubChem CID — get an inline 3D preview
and one-click downloads (FBX, OBJ+MTL, interactive HTML, settings JSON).
The resolver prefers PubChem's pre-computed 3D SDFs and falls back to
PubChem SMILES → ETKDG, then cirpy, then a direct SMILES parse.

### Run locally

```
# Pinned, reproducible:
pip install -r requirements-finder.txt
python -m streamlit run streamlit_finder.py
```

### Run in Docker

```
docker build -t molecule-finder .
docker run --rm -p 8501:8501 molecule-finder
```

The container is non-root, ships a `/_stcore/health` healthcheck, and reads
all knobs from the environment (`MIF_*` — see `finder/config.py`):

| Var                       | Default                                   |
| ------------------------- | ----------------------------------------- |
| `MIF_PUBCHEM_BASE`        | https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound |
| `MIF_HTTP_TIMEOUT`        | 15 (seconds)                              |
| `MIF_MAX_QUERY_LENGTH`    | 200                                       |
| `MIF_MAX_ATOMS`           | 600                                       |
| `MIF_MAX_BONDS`           | 600                                       |
| `MIF_CACHE_TTL`           | 86 400 (seconds)                          |
| `MIF_DEFAULT_RESOLUTION`  | 24 (FBX/OBJ triangulation)                |
| `MIF_PREVIEW_RESOLUTION`  | 20 (in-page preview)                      |
| `MIF_LOG_LEVEL`           | INFO                                      |

### Test

```
pip install -r requirements-finder.txt
pytest tests/ -v
```

30 unit + integration tests cover the mesh primitives, FBX/OBJ writers,
resolver validation, and the live PubChem branch. CI runs the same suite
on Python 3.11 and 3.12 (`.github/workflows/test-finder.yml`).

### What's inside

```
finder/
  config.py       # env-driven runtime configuration + presets
  resolver.py     # query → ResolveResult (PubChem 3D → SMILES → cirpy)
  exporters.py    # FBX / OBJ.zip / HTML / settings-JSON byte builders
mesh_export.py    # pure-Python sphere/cylinder/FBX/OBJ writers (no SDK)
streamlit_finder.py  # Streamlit UI — composition only
tests/            # offline unit tests + live integration tests
```

Quick-pick row covers THC, CBD, caffeine, aspirin, paracetamol, glucose,
adrenaline, serotonin, dopamine, nicotine, cholesterol and vitamin C.
Customization sidebar offers CPK / monochrome / fully custom palettes,
atom radius, bond thickness, mesh resolution, and the `remove hydrogens`
toggle.

## Donate

I enjoy working on this project in my free time, especially at night. If you want to support me with a coffee, just [click here!](https://www.paypal.com/donate/?hosted_button_id=V4LJ3Z3B3KXRY)

## License

Code is licensed under the GNU General Public License v3.0 ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html))

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%20v3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Graphics are licensed under the Creative Commons Non Commercial Share Alike License 4.0 ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/))

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
