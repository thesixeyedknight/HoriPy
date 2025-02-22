# **hori\_project**

**hori\_project** is a Python package designed for analyzing biomolecular structures. It offers robust tools for parsing structural files (PDB/CIF), computing molecular interactions (such as hydrogen bonds, salt bridges, π–π stacking, cation–π, and van der Waals forces), calculating interaction energies, and identifying higher-order interaction clusters.

---

## 🚀 **Features**

- 🔍 **Parsing**: Automatically detects and reads PDB/CIF file formats.
- ⚡ **Energy Computations**: Calculates hydrogen bonds, salt bridges, van der Waals forces, and other molecular interactions.
- �� **Structural Analysis**: Computes Ramachandran angles and solvent-accessible surface area (SASA).
- 🖥️ **Higher-Order Interactions**: Supports multiprocessing to analyze complex residue interactions efficiently.

---

## 📦 **Installation**

To install **hori\_project**, clone the repository and install the required dependencies:

```bash
git clone https://github.com/yourusername/hori_project.git
cd hori_project
pip install -r requirements.txt
```

---

## ⚙️ **Usage**

### Command-Line Interface (CLI)

Run structural analysis directly from the command line:

```bash
python -m hori.analysis --input examples/1a2p.cif --pH 7.0
```

### Using Within a Python Script

You can also use **hori\_project** within your own Python scripts:

```python
from hori.analysis import Hori

# Initialize the Hori class with input file and pH level
hori_instance = Hori(filename="examples/1a2p.cif", pH=7.0)

# Access computed interactions
print(hori_instance.residue_interactions)
```

---

## 🧪 **Testing**

Unit tests are included in the `tests/` directory. Run all tests using:

```bash
pytest
```

---

## 📄 **License**

This project is licensed under the [MIT License](LICENSE).

---

