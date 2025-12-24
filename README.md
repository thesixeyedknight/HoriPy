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

## Output Format

You can obtain the analysis results as a dictionary using `build_output_data` or save it directly to a JSON file using `save_output`.

```python
from hori.output import save_output, build_output_data

# Save to JSON file
save_output(hori_instance, 'output.json')

# Get dictionary directly
data = build_output_data(hori_instance)
```

The output dictionary has the following structure:

### 1. Top-Level Keys
`['meta', 'atoms', 'residues', 'atomic_interactions', 'residue_interactions', 'higher_order_interactions']`

### 2. Meta (`meta`)
Contains run metadata.
*   Keys: `['version', 'fn', 'ft', 'ho', 'pH', 'co', 'da', 'dha', 'sb', 'ppd', 'ppa', 'cpd', 'ch', 'ts']`

### 3. Atoms (`atoms`)
Keyed by integer atom IDs.
*   **Keys**: `['id', 'oid', 'resi', 'resn', 'chain', 'name', 'atom_type', 'x', 'y', 'z', 'charge', 'radius']`

### 4. Residues (`residues`)
Keyed by `{Chain}_{ResidueNumber}` (e.g., `A_24`).
*   **Keys**: `['chain', 'resi', 'resn', 'sasa_total', 'sasa_polar', 'sasa_apolar', 'sasa_main_chain', 'sasa_side_chain', 'phi', 'psi', 'omega', 'aids']`

### 5. Atomic Interactions (`atomic_interactions`)
Keyed by `{AtomID1}_{AtomID2}` (e.g., `1206_1338`).
*   **Keys**: `['a1', 'a2', 'energy', 'kbp_energy', 'type', 'dist', 'geom_metrics', 'nis']`
*   **Example**:
    ```json
    "129_1648": {
      "a1": 129,
      "a2": 1648,
      "energy": -4.2,
      "kbp_energy": -0.89,
      "type": "hbond_strong",
      "dist": 2.15,
      "geom_metrics": { ... },
      "nis": { ... }
    }
    ```

### 6. Residue Interactions (`residue_interactions`)
Keyed by `{Res1}__{Res2}` (e.g., `A_78__A_8`).
*   **Keys**: `['atomic_interaction_keys', 'int_types', 'total_energy']`
    *   `atomic_interaction_keys`: List of atomic interaction keys (strings).
    *   `int_types`: List of interaction types (e.g., `['vdw', 'salt_bridge']`).
    *   `total_energy`: Float.

### 7. Higher Order Interactions (`higher_order_interactions`)
Keyed by the order of interaction as a string (e.g., `'3'`, `'4'`).
*   Values are lists of lists, where each inner list represents an interaction group containing tuples of `(Chain, ResidueNumber)`.
*   **Example** (Quadruplet): `[[('A', 10), ('A', 88), ('A', 91), ('A', 100)], ...]`

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

