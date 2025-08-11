# hori/kbp_tools.py
import numpy as np
from pathlib import Path
import math

DATA_DIR = Path(__file__).parent.parent / "data"
KBP_FILE = DATA_DIR / "kbp_results.npz"

class KBPManager:
	"""
	Manages loading and accessing knowledge-based potentials from an NPZ file.
	"""
	def __init__(self, kbp_path=None):
		if kbp_path is None:
			kbp_path = KBP_FILE
		self.potentials = {}
		self.bin_edges = {}
		self.load_potentials(kbp_path)

	def load_potentials(self, kbp_path):
		"""
		Loads KBP data from the specified NPZ file.
		"""
		try:
			with np.load(kbp_path, allow_pickle=True) as data:
				for key in data.keys():
					if key.endswith('_edges'):
						if 'disulfide_binning_chi_edges' in key:
							self.bin_edges['disulfide_chi_ss'] = data[key]
							self.bin_edges['disulfide_chi1_cys1'] = data[key]
							self.bin_edges['disulfide_chi1_cys2'] = data[key]
						else:
							prefix = key.replace('_binning', '').replace('_edges', '')
							self.bin_edges[prefix] = data[key]
					elif "_U_kT" in key:
						prefix = key.replace('_U_kT_', '_')
						self.potentials[prefix] = data[key]
		except FileNotFoundError:
			print(f"Warning: KBP file not found at {kbp_path}. KBP energies will not be available.")
			self.potentials = {}
			self.bin_edges = {}

	def get_potential(self, interaction_type, geometry_type, value):
		"""
		Retrieves the KBP for a given interaction type, geometry, and value.
		"""
		key = f"{interaction_type}_{geometry_type}"
		potential_array = self.potentials.get(key)
		edges = self.bin_edges.get(key)
		if potential_array is None or edges is None:
			return None

		bin_index = np.digitize(value, edges) - 1

		if 0 <= bin_index < len(potential_array):
			return potential_array[bin_index]
		else:
			return None