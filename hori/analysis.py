# hori/analysis.py

# ... other imports ...
from hori.parsing import read_file, parse_pqr_atoms_list, atom_dict_to_pdb # Correctly added atom_dict_to_pdb
from hori.pdb2pqr_runner import run_pdb2pqr
from hori.interactions import (populate_bonds, build_distance_map_parallel,
							   find_atomic_interactions, find_residue_interactions)
from hori.sasa import compute_sasa_data
from hori.higher_order import parallel_find_higher_order_interactions
from hori.amber_parser import parse_amber_atomtypes, parse_amber_ffnonbonded
from hori.ramachandran import compute_ramachandran_angles

# Import for Wisz method
from hori.wisz_parameters import (WISZ_HELLINGA_PARAMS, 
							  REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED,
							  DEFAULT_IONIC_STRENGTH, DEFAULT_TEMPERATURE)

#for pike method
from hori.pike_nanda_parameters import PIKE_NANDA_ALL_PARAMS
from hori.spatial_utils import create_spatial_grid

class Hori:
	def __init__(self, filename, file_type=None, highest_order='converge',
				 pH=7.0, cutoff=7.0, d_a_dist=3.9, h_a_dist=2.5, dha_angle=90.0,
				 salt_bridge_dist=4.0, pi_pi_dist=6.0, cation_pi_dist=6.0,
				 pi_pi_angle=30.0, disulfide_min_dist=1.8, disulfide_max_dist=2.2,
				 chi_ss_angle_opts=None, chi1_angle_opts=None,
				 vdw_overlap_tolerance=0.5,
				 dielectric_method='bulk',
				 bulk_dielectric = 4.0,
				 ionic_strength=DEFAULT_IONIC_STRENGTH, 
				 temperature=DEFAULT_TEMPERATURE
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
		self.chi1_angle_opts = chi1_angle_opts if chi1_angle_opts is not None else [(-60.0, 20.0), (60.0, 20.0), (180.0, 20.0)]
		self.vdw_overlap_tolerance = vdw_overlap_tolerance

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
		self.ionic_strength = ionic_strength
		self.temperature = temperature
		self.live_ionizable_atom_sasa = {}
		self.fractional_sasa_groups = {}
		self.region_classification = {}
		self.wisz_params = None
		self.ref_sasa_ionizable = None
		if self.dielectric_method == 'wisz':
			self.wisz_params = WISZ_HELLINGA_PARAMS
			self.ref_sasa_ionizable = REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED
		
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
					# Optionally fallback, or let it try and potentially fail later if params are truly essential
			except ImportError:
				print("Error: Could not import Pike & Nanda parameters. Defaulting dielectric method to 'bulk'.")
				self.dielectric_method = 'bulk'
				self.pike_nanda_params = None

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
				grid_cell_size = influence_radius # Or max(influence_radius, some_default_like_5_or_10)
				grid, min_c, actual_cs = create_spatial_grid(self.atoms, grid_cell_size)
				self.spatial_grid_data = {'grid': grid, 'min_coords': min_c, 'cell_size': actual_cs}
				print(f"Spatial grid created for Pike & Nanda method with cell size: {actual_cs:.2f} Å")
			else:
				print("Warning: POLARIZABILITY_INFLUENCE_RADIUS not defined or invalid in Pike & Nanda params. Spatial grid not created. Pike method may fail.")
				self.spatial_grid_data = None # Ensure it's None
				# Consider falling back to bulk if grid is essential and cannot be created
				# self.dielectric_method = 'bulk' 

		# --- SASA Calculation ---
		sasa_calculation_successful = False
		all_atom_sasa_map_from_fs = {}
		if self.atoms:
			# This import is fine here, or at the top of the file.
			from hori.wisz_utils import _get_ionizable_group_type_for_sasa 
			
			pdb_lines_for_sasa_str = atom_dict_to_pdb(self.atoms)
			input_lines_for_sasa_calc = [line + "\n" for line in pdb_lines_for_sasa_str.split('\n') if line.strip()]

			sasa_error, all_atom_sasa_map_from_fs = compute_sasa_data(
				input_lines_for_sasa_calc,
				self.residues, 
				'pdb' 
			)

			if sasa_error:
				print("Error during SASA calculation. Subsequent Wisz calculations will be affected.")
			else:
				sasa_calculation_successful = True
				if self.dielectric_method == 'wisz' and self.ref_sasa_ionizable: # Check ref_sasa_ionizable
					for (atom_fs_chain, atom_fs_resi, atom_fs_name), sasa_val in all_atom_sasa_map_from_fs.items():
						res_record = self.residues.get((atom_fs_chain, atom_fs_resi))
						if res_record:
							hori_resn = res_record.resn.upper()
							group_type = _get_ionizable_group_type_for_sasa(
								hori_resn, atom_fs_name, self.ref_sasa_ionizable
							)
							if group_type:
								self.live_ionizable_atom_sasa[(atom_fs_chain, atom_fs_resi, atom_fs_name)] = sasa_val
		else:
			print("Warning: No Hori atoms available for SASA calculation.")
		if self.dielectric_method == 'wisz':
			if sasa_calculation_successful and self.live_ionizable_atom_sasa:
				from hori.wisz_utils import calculate_fractional_sasa_for_groups
				self.fractional_sasa_groups = calculate_fractional_sasa_for_groups(
					self.live_ionizable_atom_sasa,
					self.residues,
					self.ref_sasa_ionizable
				)
			elif not sasa_calculation_successful:
				print("Warning: SASA calculation failed. Fractional SASA for Wisz method cannot be calculated.")
			else: 
				print("Warning: Live ionizable atom SASA dictionary is empty (no Wisz-relevant ionizable atoms found/mapped). "
						"Fractional SASA for Wisz method cannot be calculated.")
		# --- End SASA Calculation ---
		
		# --- Region Classification (Only if Wisz method is active and prerequisites met) ---
		if self.dielectric_method == 'wisz' and sasa_calculation_successful and self.fractional_sasa_groups:
			from hori.wisz_utils import calculate_cb_depths, classify_residue_region
			from hori.parsing import prepare_pdb_for_pseudo_surface
			
			pdb_lines_bb_cb, bb_cb_coords_list = prepare_pdb_for_pseudo_surface(self.atoms)
			cb_depths = {} 
			if pdb_lines_bb_cb and bb_cb_coords_list: 
				cb_depths = calculate_cb_depths(
							self.atoms,      
							self.residues,   
							pdb_lines_bb_cb, 
							bb_cb_coords_list
							)
			else:
				print("Warning: Skipping C-beta depth calculation as no backbone/CB atoms were prepared for pseudo-surface.")
				for (chain, resi) in self.residues.keys():
					cb_depths[(chain, resi)] = 99.0 

			self.region_classification = {} 
			for group_key, af_value in self.fractional_sasa_groups.items():
				chain, resi, group_type_str = group_key 
				r_value = cb_depths.get((chain, resi))

				if r_value is not None:
					region = classify_residue_region(r_value, af_value)
					self.region_classification[group_key] = region
				else:
					print(f"Warning: C-beta depth for residue {chain}{resi} not found. Region for group {group_type_str} set to 'surface'.")
					self.region_classification[group_key] = 'surface' 
		elif self.dielectric_method == 'wisz':
			print("Warning: Skipping Wisz region classification due to missing SASA data or fractional SASA groups.")
		# --- End Region Classification ---

		# --- Setup user_params (should be done after all Wisz-specific data is computed) ---
		self.user_params = {
			'd_a_dist': self.d_a_dist, 'dha_angle': self.dha_angle, 'h_a_dist': self.h_a_dist,
			'salt_bridge_dist': self.salt_bridge_dist, 'pi_pi_dist': self.pi_pi_dist,
			'cation_pi_dist': self.cation_pi_dist, 'pi_pi_angle': self.pi_pi_angle,
			'disulfide_min_dist': self.disulfide_min_dist, 'disulfide_max_dist': self.disulfide_max_dist,
			'chi_ss_angle_opts': self.chi_ss_angle_opts, 'chi1_angle_opts': self.chi1_angle_opts,
			'vdw_overlap_tolerance': self.vdw_overlap_tolerance,
			'dielectric_method': self.dielectric_method
		}

		if self.dielectric_method == 'bulk':
			self.user_params['bulk_dielectric_value'] = self.bulk_dielectric_value
		elif self.dielectric_method == 'wisz':
			if self.wisz_params and self.ref_sasa_ionizable:
				self.user_params['wisz_params'] = self.wisz_params
				self.user_params['fractional_sasa_groups'] = self.fractional_sasa_groups # Now populated
				self.user_params['region_classification'] = self.region_classification # Now populated
				self.user_params['ionic_strength'] = self.ionic_strength
				self.user_params['temperature'] = self.temperature
				self.user_params['ref_sasa_ionizable'] = self.ref_sasa_ionizable
			else:
				print("Warning: Wisz method selected, but Wisz parameters are not fully initialized. Interaction energies may be incorrect.")
		elif self.dielectric_method == 'pike_nanda': # (NEW Block)
			if self.pike_nanda_params and self.spatial_grid_data: # Check if both are successfully initialized
				self.user_params['pike_nanda_params'] = self.pike_nanda_params
				self.user_params['spatial_grid_data'] = self.spatial_grid_data
			else:
				print(f"Warning: Pike & Nanda method selected, but parameters or spatial grid are not initialized. Falling back to bulk dielectric for energy calculations.")
				self.user_params['dielectric_method'] = 'bulk' # Fallback for energy calculations
				self.user_params['bulk_dielectric_value'] = self.bulk_dielectric_value # Ensure bulk value is available

# --- End Setup user_params ---

		# --- Interaction Calculations ---
		build_distance_map_parallel(self.residues, self.distance_map, cutoff=self.cutoff)
		find_atomic_interactions(self.distance_map, self.atoms, self.bonds, self.residues,
								 self.atom_interactions, self.amber_nonbonded, self.user_params, self)
		find_residue_interactions(self.atom_interactions, self.residue_interactions)

		self.higher_order_interactions = parallel_find_higher_order_interactions(
			self.residues, self.residue_interactions, max_order=self.highest_order
		)