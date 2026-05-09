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

## Donate

I enjoy working on this project in my free time, especially at night. If you want to support me with a coffee, just [click here!](https://www.paypal.com/donate/?hosted_button_id=V4LJ3Z3B3KXRY)

## License

Code is licensed under the GNU General Public License v3.0 ([GPL-3.0](https://www.gnu.org/licenses/gpl-3.0.en.html))

[![License: GPL-3.0](https://img.shields.io/badge/License-GPL%20v3-lightgrey.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

Graphics are licensed under the Creative Commons Non Commercial Share Alike License 4.0 ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/))

[![License: CC BY-NC-SA 4.0](https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by-nc-sa/4.0/)
