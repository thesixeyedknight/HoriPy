# Hori: Molecular Interaction Analysis Toolkit

**Hori** is a Python library for analyzing molecular structures, focusing on computing higher-order interactions, solvent-accessible surface area (SASA), and energy-based residue interactions in biomolecules.

## Features
- Parses **PDB** and **CIF** files for protein structures.
- Computes atomic and residue-level interactions (hydrogen bonds, salt bridges, van der Waals forces, etc.).
- Calculates solvent-accessible surface area (SASA) using `freesasa`.
- Generates higher-order residue interaction networks.
- Supports customizable interaction parameters like pH, distance cutoffs, and angles.

## Installation

Clone the repository and install dependencies:

```bash
git clone https://github.com/thesixeyedknight/hori.git
cd hori
pip install -r requirements.txt
```

## Usage Example

```python
import sys
sys.path.append('/path/to/hori_project')

from hori.analysis import Hori
from hori.output import save_output

# Load a structure and perform analysis
h1 = Hori(filename='/path/to/1BVP.cif', pH=7.0)
save_output(h1, '1bvp_output.json')

# Access computed residue data
for key, residue in h1.residues.items():
    print(f"Residue {key}: SASA Total = {residue.sasa_total:.2f}")
```

## Requirements

- `propka`
- `numpy`
- `freesasa`
- `tqdm`
- `requests`

Install them with:

```bash
pip install -r requirements.txt
```

## Running Tests

```bash
python -m unittest discover tests
```

## License

- The majority of this repository is licensed under the GNU GPLv3. See the [LICENSE](./LICENSE) file for details.
- The code in the `mod_pdb2pqr` directory is based on [pdb2pqr](https://github.com/Electrostatics/pdb2pqr), which is distributed under its original license. See the [license file in that directory](./mod_pdb2pqr/LICENCE) for more information.

