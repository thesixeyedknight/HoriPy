# hori/analysis.py

# ... method imports ...
from hori.parsing import read_file, parse_pqr_atoms_list, atom_dict_to_pdb
from hori.pdb2pqr_runner import run_pdb2pqr
from hori.interactions import (populate_bonds, build_distance_map_parallel,
							   find_atomic_interactions, find_residue_interactions)
from hori.sasa import compute_sasa_data
from hori.higher_order import parallel_find_higher_order_interactions
from hori.amber_parser import parse_amber_atomtypes, parse_amber_ffnonbonded
from hori.ramachandran import compute_ramachandran_angles
from hori.kbp_tools import KBPManager # for knowledge based potentials
#for pike method
from hori.pike_nanda_parameters import PIKE_NANDA_ALL_PARAMS
from hori.spatial_utils import create_spatial_grid

class Hori:
	def __init__(self, filename, file_type=None, highest_order='converge',
				 pH=7.0, cutoff=7.0, d_a_dist=3.9, h_a_dist=2.5, dha_angle=90.0,
				 salt_bridge_dist=4.0, pi_pi_dist=6.0, cation_pi_dist=6.0,
				 pi_pi_angle=30.0, disulfide_min_dist=1.7, disulfide_max_dist=2.4,
				 chi_ss_angle_opts=None, chi1_angle_opts=None,
				 vdw_overlap_tolerance=0.5,
				 dielectric_method='bulk',
				 bulk_dielectric = 4.0,
				 num_cores=8,
				 use_kbp=True,
				 kbp_file_path=None
				 ):

		# ... (parameter storage) ...
		self.filename = filename
		self.chains = set()
		self.highest_order = highest_order
		self.ph = float(pH)
		self.cutoff = cutoff
		self.d_a_dist = d_a_dist
		self.h_a_dist = h_a_dist
		self.dha_angle = dha_angle
		self.salt_bridge_dist = salt_bridge_dist
		self.pi_pi_dist = pi_pi_dist
		self.pi_pi_angle = pi_pi_angle
		self.cation_pi_dist = cation_pi_dist
		self.disulfide_min_dist = disulfide_min_dist
		self.disulfide_max_dist = disulfide_max_dist
		self.chi_ss_angle_opts = chi_ss_angle_opts if chi_ss_angle_opts is not None else [(97.0, 30.0), (-87.0, 30.0)]
		self.chi1_angle_opts = chi1_angle_opts if chi1_angle_opts is not None else [(-60.0, 35.0), (60.0, 35.0), (180.0, 35.0)]
		self.vdw_overlap_tolerance = vdw_overlap_tolerance
		self.num_cores = num_cores

		self.atoms = {}
		self.residues = {}
		self.bonds = {}
		self.pdb_for_haad = []
		self.atom_interactions = {}
		self.distance_map = {}
		self.amber_atomtypes = {}
		self.amber_nonbonded = {}
		self.residue_interactions = {}
		self.pqr_atoms = []

		self.dielectric_method = dielectric_method
		self.bulk_dielectric_value = float(bulk_dielectric)
		
		self.pike_nanda_params = None
		self.spatial_grid_data = None
		if self.dielectric_method == 'pike_nanda':
			try:
				self.pike_nanda_params = PIKE_NANDA_ALL_PARAMS
				if not self.pike_nanda_params.get('PIKE_ATOMIC_POLARIZABILITIES') or \
				   not self.pike_nanda_params.get('PIKE_ATOMIC_VOLUMES') or \
				   not self.pike_nanda_params.get('AMBER_TO_PIKE_ATOM_TYPE_MAP') or \
				   self.pike_nanda_params.get('POLARIZABILITY_INFLUENCE_RADIUS') is None:
					print("Warning: Core Pike & Nanda parameters (polarizabilities, volumes, map, radius) missing. Pike & Nanda method may fail or be inaccurate.")
			except ImportError: # Should be caught by the initial import if file is missing
				print("Error: Could not import Pike & Nanda parameters. Defaulting dielectric method to 'bulk'.")
				self.dielectric_method = 'bulk'
				self.pike_nanda_params = None

		self.kbp_manager = None
		if use_kbp:
			self.kbp_manager = KBPManager(kbp_file_path)

		# --- Workflow ---
		parse_amber_atomtypes(self.amber_atomtypes)
		parse_amber_ffnonbonded(self.amber_nonbonded)

		self.file_type, self.pdb_for_haad, self.chains = read_file(self.filename)
		print(f'File format: {self.file_type}')

		self.pqr_atoms = run_pdb2pqr(self.pdb_for_haad, self.file_type, self.ph)
		if not self.pqr_atoms:
			print("Error: PDB2PQR execution failed or returned no atoms. Aborting Hori analysis.")
			return

		self.atoms, self.residues, self.bonds = parse_pqr_atoms_list(self.pqr_atoms)
		populate_bonds(self.residues, self.bonds)
		self.residues = compute_ramachandran_angles(self.residues)

		# --- pike nanda calculations ---
		# --- Create Spatial Grid (if using Pike & Nanda) ---
		if self.dielectric_method == 'pike_nanda' and self.pike_nanda_params and self.atoms:
			influence_radius = self.pike_nanda_params.get('POLARIZABILITY_INFLUENCE_RADIUS')
			if influence_radius is not None and influence_radius > 0:
				grid_cell_size = influence_radius 
				grid, min_c, actual_cs = create_spatial_grid(self.atoms, grid_cell_size)
				self.spatial_grid_data = {'grid': grid, 'min_coords': min_c, 'cell_size': actual_cs}
				print(f"Spatial grid created for Pike & Nanda method with cell size: {actual_cs:.2f} Å")
			else:
				print("Warning: POLARIZABILITY_INFLUENCE_RADIUS not defined or invalid in Pike & Nanda params. Spatial grid not created. Pike method may fail.")
				self.spatial_grid_data = None 
				# Fallback to bulk if grid is essential
				print("Falling back to 'bulk' dielectric method as spatial grid for Pike & Nanda could not be created.")
				self.dielectric_method = 'bulk'


		# --- SASA Calculation ---
		# sasa_calculation_successful = False # Not needed if not used for Wisz
		all_atom_sasa_map_from_fs = {}
		if self.atoms:
			pdb_lines_for_sasa_str = atom_dict_to_pdb(self.atoms)
			input_lines_for_sasa_calc = [line + "\n" for line in pdb_lines_for_sasa_str.split('\n') if line.strip()]

			sasa_error, all_atom_sasa_map_from_fs = compute_sasa_data(
				input_lines_for_sasa_calc,
				self.residues,
				'pdb'
			)
			if sasa_error:
				print("Error during SASA calculation.")
			# else:
				# sasa_calculation_successful = True # Not explicitly needed without Wisz logic following immediately
		else:
			print("Warning: No Hori atoms available for SASA calculation.")
		# --- End SASA Calculation ---

		# --- Setup user_params ---
		self.user_params = {
			'd_a_dist': self.d_a_dist, 'dha_angle': self.dha_angle, 'h_a_dist': self.h_a_dist,
			'salt_bridge_dist': self.salt_bridge_dist, 'pi_pi_dist': self.pi_pi_dist,
			'cation_pi_dist': self.cation_pi_dist, 'pi_pi_angle': self.pi_pi_angle,
			'disulfide_min_dist': self.disulfide_min_dist, 'disulfide_max_dist': self.disulfide_max_dist,
			'chi_ss_angle_opts': self.chi_ss_angle_opts, 'chi1_angle_opts': self.chi1_angle_opts,
			'vdw_overlap_tolerance': self.vdw_overlap_tolerance,
			'dielectric_method': self.dielectric_method # This will be 'bulk' if pike_nanda failed above
		}

		if self.dielectric_method == 'bulk':
			self.user_params['bulk_dielectric_value'] = self.bulk_dielectric_value
		elif self.dielectric_method == 'pike_nanda':
			# This check is important: if spatial_grid_data wasn't created,
			# dielectric_method would have been reset to 'bulk'.
			if self.pike_nanda_params and self.spatial_grid_data:
				self.user_params['pike_nanda_params'] = self.pike_nanda_params
				self.user_params['spatial_grid_data'] = self.spatial_grid_data
			else:
				# This case should ideally be prevented by the fallback logic above,
				# but as a safeguard:
				print("Error: Pike & Nanda selected but prerequisites missing. Reverting to bulk dielectric for safety.")
				self.user_params['dielectric_method'] = 'bulk'
				self.user_params['bulk_dielectric_value'] = self.bulk_dielectric_value


		# --- Interaction Calculations ---
		build_distance_map_parallel(self.residues, self.distance_map, cutoff=self.cutoff, num_cores=self.num_cores)
		find_atomic_interactions(self.distance_map, self.atoms, self.bonds, self.residues,
								 self.atom_interactions, self.amber_nonbonded, self.user_params, self, kbp_manager=self.kbp_manager)
		find_residue_interactions(self.atom_interactions, self.residue_interactions)

		self.higher_order_interactions = parallel_find_higher_order_interactions(
			self.residues, self.residue_interactions, max_order=self.highest_order
		)