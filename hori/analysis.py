from hori.parsing import read_file, parse_pqr_atoms_list #, read_pdb_with_hydrogen
from hori.pdb2pqr_runner import run_pdb2pqr #, run_pdb2pqr_old
from hori.interactions import populate_bonds, build_distance_map_parallel, find_atomic_interactions, find_residue_interactions
from hori.sasa import compute_residue_sasa
from hori.higher_order import parallel_find_higher_order_interactions
from hori.amber_parser import parse_amber_atomtypes, parse_amber_ffnonbonded
from hori.ramachandran import compute_ramachandran_angles

class Hori:

	def __init__(self, filename, file_type=None, highest_order='converge',
			  pH=7.0, cutoff=7.0, d_a_dist=3.2, dha_angle=150.0, salt_bridge_dist=4.0,
			  pi_pi_dist=6.0, cation_pi_dist=6.0, pi_pi_angle=30.0):
		self.filename = filename
		self.file_type = file_type
		self.chains = set()
		self.highest_order = highest_order
		self.ph = float(pH)
		self.cutoff = cutoff
		self.d_a_dist = d_a_dist
		self.dha_angle = dha_angle
		self.salt_bridge_dist = salt_bridge_dist
		self.pi_pi_dist = pi_pi_dist
		self.pi_pi_angle = pi_pi_angle
		self.cation_pi_dist = cation_pi_dist
		self.atoms = {}
		self.residues = {}
		self.bonds = {}
		self.pdb_for_haad = []
		self.pdb_with_hydrogen = []
		self.atom_interactions = {}
		self.distance_map = {}
		self.amber_atomtypes = {}
		self.amber_nonbonded = {}
		self.residue_interactions = {}

		# Parse AMBER data
		parse_amber_atomtypes(self.amber_atomtypes)
		parse_amber_ffnonbonded(self.amber_nonbonded)

		# Read input file
		self.file_type, self.pdb_for_haad, self.chains = read_file(self.filename)
		print(f'File format: {self.file_type}')
		# Run pdb2pqr
		self.pqr_atoms = run_pdb2pqr(self.pdb_for_haad, self.file_type, self.ph)
		#self.pdb_with_hydrogen = run_pdb2pqr_old(self.pdb_for_haad)
		
		#read pdb with hydrogens
		#self.atoms, self.residues, self.bonds = read_pdb_with_hydrogen(self.pdb_with_hydrogen)
		self.atoms, self.residues, self.bonds = parse_pqr_atoms_list(self.pqr_atoms)

		#build bond map
		populate_bonds(self.residues, self.bonds)

		# Compute Ramachandran angles
		self.residues = compute_ramachandran_angles(self.residues)

		# Compute SASA
		compute_residue_sasa(self.pdb_for_haad, self.residues, self.file_type)

		# Compute interactions
		build_distance_map_parallel(self.residues, self.distance_map, cutoff=self.cutoff)

		self.user_params = {'d_a_dist': self.d_a_dist, 'dha_angle': self.dha_angle,
					   'salt_bridge_dist': self.salt_bridge_dist, 'pi_pi_dist': self.pi_pi_dist,
					   'cation_pi_dist': self.cation_pi_dist, 'pi_pi_angle': self.pi_pi_angle}
		
		find_atomic_interactions(self.distance_map, self.atoms, self.bonds, self.residues, self.atom_interactions, self.amber_nonbonded, self.user_params)
		find_residue_interactions(self.atom_interactions, self.residue_interactions)

		# Higher-order interactions
		self.higher_order_interactions = parallel_find_higher_order_interactions(self.residues, self.residue_interactions, max_order=self.highest_order)
