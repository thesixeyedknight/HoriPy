import unittest
import os
from hori.analysis import Hori

class TestHori(unittest.TestCase):

	def setUp(self):
		# Paths to test data and required files
		self.example_pdb = os.path.join("examples", "1a2p.cif")
		#self.atomtypes_file = os.path.join("data", "atomtypes.atp")
		#self.nonbonded_file = os.path.join("data", "ffnonbonded.itp")

		# Ensure the files exist
		for file in [self.example_pdb]:
			self.assertTrue(os.path.exists(file), f"Required file missing: {file}")

	def test_hori_initialization(self):
		"""Test if the Hori class initializes correctly."""
		hori_instance = Hori(filename=self.example_pdb)
		self.assertIsNotNone(hori_instance)
		self.assertEqual(hori_instance.filename, self.example_pdb)
		self.assertTrue(hori_instance.atoms)
		self.assertTrue(hori_instance.residues)

	def test_sasa_computation(self):
		"""Test if SASA is computed and added to residues."""
		hori_instance = Hori(filename=self.example_pdb)
		for residue in hori_instance.residues.values():
			self.assertGreaterEqual(residue.sasa_total, 0)

	def test_interactions(self):
		"""Test interaction computations."""
		hori_instance = Hori(filename=self.example_pdb)
		self.assertTrue(hori_instance.atom_interactions)
		self.assertTrue(hori_instance.residue_interactions)

	def test_hori(self):
		hori_instance = Hori(filename=self.example_pdb)
		self.assertTrue(hori_instance.higher_order_interactions)

if __name__ == "__main__":
	unittest.main()
