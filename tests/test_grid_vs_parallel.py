"""
Comparison test: grid-based distance map vs. original parallel distance map.

Runs both approaches on the same PDB and asserts they produce identical results.
"""
import unittest
import os
import math

from hori.parsing import read_file, parse_pqr_atoms_list
from hori.pdb2pqr_runner import run_pdb2pqr
from hori.interactions import (
	populate_bonds, build_distance_map_parallel,
	build_distance_map_from_grid
)
from hori.spatial_utils import create_spatial_grid
from hori.ramachandran import compute_ramachandran_angles


class TestGridVsParallel(unittest.TestCase):

	@classmethod
	def setUpClass(cls):
		"""Parse the test PDB once for all tests."""
		example_pdb = os.path.join("examples", "1a2p.cif")
		if not os.path.exists(example_pdb):
			raise FileNotFoundError(f"Test file missing: {example_pdb}")

		file_type, pdb_for_haad, chains = read_file(example_pdb)
		pqr_atoms = run_pdb2pqr(pdb_for_haad, file_type, 7.0)
		cls.atoms, cls.residues, cls.bonds = parse_pqr_atoms_list(pqr_atoms)
		populate_bonds(cls.residues, cls.bonds)
		cls.residues = compute_ramachandran_angles(cls.residues)

		cls.cutoff = 7.0

		# Build distance map using the ORIGINAL parallel method
		cls.parallel_distance_map = {}
		build_distance_map_parallel(
			cls.residues, cls.parallel_distance_map,
			cutoff=cls.cutoff, num_cores=1
		)

		# Build distance map using the NEW grid-based method
		grid, min_c, actual_cs = create_spatial_grid(cls.atoms, cls.cutoff)
		grid_data = {'grid': grid, 'min_coords': min_c, 'cell_size': actual_cs}
		cls.grid_distance_map = build_distance_map_from_grid(
			cls.atoms, cls.residues, grid_data, cutoff=cls.cutoff
		)

	def test_same_keys(self):
		"""Both methods should find the same atom pairs."""
		parallel_keys = set(self.parallel_distance_map.keys())
		grid_keys = set(self.grid_distance_map.keys())

		# Keys only in parallel but not grid
		only_parallel = parallel_keys - grid_keys
		# Keys only in grid but not parallel
		only_grid = grid_keys - parallel_keys

		self.assertEqual(
			only_parallel, set(),
			f"{len(only_parallel)} pairs found by parallel but not grid"
		)
		self.assertEqual(
			only_grid, set(),
			f"{len(only_grid)} pairs found by grid but not parallel"
		)

	def test_same_distances(self):
		"""Distances should match within floating-point tolerance."""
		tolerance = 1e-9
		mismatched = []
		for key in self.parallel_distance_map:
			if key in self.grid_distance_map:
				d_par = self.parallel_distance_map[key]
				d_grid = self.grid_distance_map[key]
				if abs(d_par - d_grid) > tolerance:
					mismatched.append((key, d_par, d_grid))

		self.assertEqual(
			len(mismatched), 0,
			f"{len(mismatched)} pairs have distance mismatches (tolerance={tolerance})"
		)

	def test_map_sizes_equal(self):
		"""Both maps should have the same number of entries."""
		self.assertEqual(
			len(self.parallel_distance_map),
			len(self.grid_distance_map),
			f"Size mismatch: parallel={len(self.parallel_distance_map)}, "
			f"grid={len(self.grid_distance_map)}"
		)


if __name__ == "__main__":
	unittest.main()
