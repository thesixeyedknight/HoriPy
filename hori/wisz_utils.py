# hori/wisz_utils.py
import math
import tempfile
import freesasa as fs
# Assuming wisz_parameters.py is in the same directory (hori/)
from hori.wisz_parameters import REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED # Already there
from hori.wisz_parameters import DEFAULT_WATER_DIELECTRIC, DEFAULT_IONIC_STRENGTH, DEFAULT_TEMPERATURE

fs.setVerbosity(1)

def _get_ionizable_group_type_for_sasa(residue_name, atom_name, ref_sasa_ionizable_dict):
	"""
	Helper function to determine the ionizable group type (e.g., "ASP_SIDECHAIN")
	for a given atom based on the reference SASA dictionary.
	This is specifically for grouping atoms correctly for fractional SASA calculation.

	Args:
		residue_name (str): The 3-letter code of the residue (e.g., "ASP").
		atom_name (str): The name of the atom (e.g., "OD1").
		ref_sasa_ionizable_dict (dict): The REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED dictionary.

	Returns:
		str or None: The group type key (e.g., "ASP_SIDECHAIN") or None if not found.
	"""
	# Handle specific cases like NTERM for PRO vs other residues
	if residue_name == "PRO" and atom_name == "N":
		if "NTERM_PRO" in ref_sasa_ionizable_dict and atom_name in ref_sasa_ionizable_dict["NTERM_PRO"]["atoms"]:
			return "NTERM_PRO"
	
	# General cases
	for group_key, info in ref_sasa_ionizable_dict.items():
		# Extract the base residue name or terminal type from the group_key
		# e.g., "ASP_SIDECHAIN" -> "ASP", "NTERM" -> "NTERM", "NTERM_PRO" -> "NTERM_PRO"
		key_parts = group_key.split('_')
		res_part_of_key_from_dict = key_parts[0]

		# Check if the atom is part of this group
		if atom_name in info["atoms"]:
			if res_part_of_key_from_dict == residue_name: # Exact match for standard residues
				return group_key
			# Handle NTERM/CTERM cases
			elif res_part_of_key_from_dict == "NTERM":
				if group_key == "NTERM_PRO" and residue_name == "PRO":
					return "NTERM_PRO"
				elif group_key == "NTERM" and residue_name != "PRO": # General N-terminus
					return "NTERM"
			elif res_part_of_key_from_dict == "CTERM" and group_key == "CTERM": # General C-terminus
				return "CTERM"
	return None


def calculate_fractional_sasa_for_groups(live_ionizable_atom_sasa_dict,
										 hori_residues_dict,
										 reference_sasa_groups_dict):
	"""
	Calculates the fractional SASA for predefined ionizable groups in residues.

	Args:
		live_ionizable_atom_sasa_dict (dict): SASA for individual ionizable atoms.
											  Format: {(chain, resi_num, atom_name): sasa_value}
		hori_residues_dict (dict): Hori's main residues dictionary.
								   Format: {(chain, resi_num): ResidueRecord}
		reference_sasa_groups_dict (dict): REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED.

	Returns:
		dict: Fractional SASA for each group in each relevant residue.
			  Format: {(chain, resi_num, group_type_str): fractional_sasa_value}
	"""
	fractional_sasa_results = {}
	aggregated_live_sasa_for_groups = {}  # Stores sum of live SASA for atoms in a group

	# Step 1: Aggregate live SASA for atoms within each defined ionizable group per residue
	for (atom_chain, atom_resi_num, atom_name), live_sasa in live_ionizable_atom_sasa_dict.items():
		res_record = hori_residues_dict.get((atom_chain, atom_resi_num))
		if not res_record:
			# This should ideally not happen if live_ionizable_atom_sasa_dict is built correctly
			print(f"Warning: Residue record for {atom_chain}{atom_resi_num} not found while calculating fractional SASA.")
			continue
		
		residue_name_from_record = res_record.resn.upper() # Ensure uppercase for matching

		group_type = _get_ionizable_group_type_for_sasa(residue_name_from_record, atom_name, reference_sasa_groups_dict)

		if group_type:
			group_instance_key = (atom_chain, atom_resi_num, group_type)
			aggregated_live_sasa_for_groups[group_instance_key] = \
				aggregated_live_sasa_for_groups.get(group_instance_key, 0.0) + live_sasa

	# Step 2: Calculate fractional SASA using the aggregated live SASA and reference SASA
	for (chain, resi_num, group_type), total_live_sasa_for_group in aggregated_live_sasa_for_groups.items():
		reference_group_info = reference_sasa_groups_dict.get(group_type)

		if reference_group_info and reference_group_info["sasa"] > 1e-9: # Avoid division by zero/small number
			frac_sasa = total_live_sasa_for_group / reference_group_info["sasa"]
			fractional_sasa_results[(chain, resi_num, group_type)] = min(frac_sasa, 1.0) # Cap at 1.0 as per typical definitions
		else:
			fractional_sasa_results[(chain, resi_num, group_type)] = 0.0 # Default if no ref or ref is zero
			if not reference_group_info:
				print(f"Warning: No reference SASA found for group type '{group_type}' "
					  f"in residue {chain}{resi_num} during fractional SASA calculation.")
			elif reference_group_info["sasa"] <= 1e-9:
				print(f"Warning: Reference SASA for group type '{group_type}' "
					  f"in residue {chain}{resi_num} is zero or very small. Fractional SASA set to 0.")
					  
	return fractional_sasa_results


def calculate_cb_depths(hori_atoms_dict_full_structure, 
						  hori_residues_dict, 
						  pdb_lines_bb_cb_only,
						  bb_cb_atom_coords_ordered_list): # New argument
	"""
	Calculates the depth (r) for the C-beta atom of each residue.
	Depth is the minimum distance from the C-beta to a pseudo-surface
	generated by a 7A probe on backbone and C-beta atoms only.

	Args:
		hori_atoms_dict_full_structure (dict): Hori's main atom dictionary (id -> Atom) for the full structure.
		hori_residues_dict (dict): Hori's main residue dictionary ((chain, resi) -> ResidueRecord).
		pdb_lines_bb_cb_only (list): List of PDB strings for backbone/CB atoms only.
		bb_cb_atom_coords_ordered_list (list): List of (x,y,z) tuples for the atoms in 
											   pdb_lines_bb_cb_only, in the same order.

	Returns:
		dict: Mapping (chain, resi) to C-beta depth (r).
	"""
	cb_depths_dict = {}
	pseudo_surface_atom_coords = []

	if not pdb_lines_bb_cb_only or not bb_cb_atom_coords_ordered_list:
		print("Warning: Cannot calculate CB depths, PDB lines or coordinates for pseudo-surface are empty.")
		for (chain, resi), _ in hori_residues_dict.items():
			cb_depths_dict[(chain, resi)] = 99.0 
		return cb_depths_dict

	pdb_text_bb_cb = "".join(pdb_lines_bb_cb_only)

	with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb') as tmp_bb_cb:
		tmp_bb_cb.write(pdb_text_bb_cb)
		tmp_bb_cb.flush()

		try:
			structure_bb_cb = fs.Structure(tmp_bb_cb.name)
			params = fs.Parameters()
			params.setProbeRadius(7.0)
			result_bb_cb = fs.calc(structure_bb_cb, params)
			#atom_sasa_bb_cb = result_bb_cb.atomAreas()

			# Check if number of atoms in FreeSASA structure matches our coordinate list
			if structure_bb_cb.nAtoms() != len(bb_cb_atom_coords_ordered_list):
				print(f"Warning: Mismatch in atom count for pseudo-surface calculation. FS: {structure_bb_cb.nAtoms()}, Hori: {len(bb_cb_atom_coords_ordered_list)}")
				# Fallback if counts don't match
				#for (chain, resi), _ in hori_residues_dict.items():
				#	cb_depths_dict[(chain, resi)] = 99.0
				#return cb_depths_dict

			for i in range(structure_bb_cb.nAtoms()):
				atom_sasa_value = result_bb_cb.atomArea(i)
				if atom_sasa_value > 0.01:
					if i < len(bb_cb_atom_coords_ordered_list):
						pseudo_surface_atom_coords.append(bb_cb_atom_coords_ordered_list[i])
					else:
						print(f"Warning: FreeSASA atom index {i} for pseudo-surface out of bounds for provided coordinate list (len {len(bb_cb_atom_coords_ordered_list)}). Skipping this surface atom.")
		
		except RuntimeError as e: # Catch FreeSASA specific runtime errors too
			print(f"Error running FreeSASA for pseudo-surface: {e}")
			for (chain, resi), _ in hori_residues_dict.items():
				cb_depths_dict[(chain, resi)] = 99.0
			return cb_depths_dict
		except Exception as e: # Generic catch
			print(f"Unexpected error during pseudo-surface SASA calculation: {e}")
			for (chain, resi), _ in hori_residues_dict.items():
				cb_depths_dict[(chain, resi)] = 99.0
			return cb_depths_dict


	if not pseudo_surface_atom_coords:
		print("Warning: No atoms found on the 7A-probe pseudo-surface. All C-beta depths will be set to 99.0.")
		for (chain, resi), _ in hori_residues_dict.items():
			cb_depths_dict[(chain, resi)] = 99.0
		return cb_depths_dict

	# Calculate depth (r) for each residue's C-beta in the original structure
	for (chain, resi_num), res_record in hori_residues_dict.items():
		# Use the hori_atoms_dict_full_structure to get Cbeta/CA coordinates
		target_atom_for_depth = None
		if res_record.resn.upper() == "GLY":
			target_atom_for_depth = next((atom for atom in res_record.atoms if atom.name.strip().upper() == "CA"), None)
		else:
			target_atom_for_depth = next((atom for atom in res_record.atoms if atom.name.strip().upper() == "CB"), None)

		if target_atom_for_depth:
			min_dist_to_surface = 999.0
			target_atom_coords = (target_atom_for_depth.x, target_atom_for_depth.y, target_atom_for_depth.z)
			for surface_coord in pseudo_surface_atom_coords:
				dist = math.sqrt(
					(target_atom_coords[0] - surface_coord[0])**2 +
					(target_atom_coords[1] - surface_coord[1])**2 +
					(target_atom_coords[2] - surface_coord[2])**2
				)
				if dist < min_dist_to_surface:
					min_dist_to_surface = dist
			cb_depths_dict[(chain, resi_num)] = min_dist_to_surface
		else:
			# print(f"Warning: CB/CA atom not found for residue {chain}{resi_num} ({res_record.resn}). Cannot calculate depth.")
			cb_depths_dict[(chain, resi_num)] = 99.0 # Default if no CB/CA
			
	return cb_depths_dict


def classify_residue_region(cb_depth_r, fractional_sasa_Af_group):
	"""
	Classifies a residue's ionizable group into Core, Boundary, or Surface
	based on its C-beta depth (r) and the group's fractional SASA (Af).

	Args:
		cb_depth_r (float): Depth of the C-beta atom from the pseudo-surface.
		fractional_sasa_Af_group (float): Fractional SASA of the ionizable group.

	Returns:
		str: 'core', 'boundary', or 'surface'.
	"""
	r = cb_depth_r
	Af = fractional_sasa_Af_group

	# Ensure Af is within [0,1] range, cap if slightly outside due to precision
	Af = max(0.0, min(Af, 1.0))


	if (r >= 7.0 and Af < 0.30) or \
	   (5.7 < r < 7.0 and Af <= 0.10): # Note: paper uses < for r < 7.0 but also r >= 7.0. Assuming 5.7 < r < 7.0 here.
									  # The criteria are r>=7 AND Af < 0.30 OR (5.7 < r < 7.0 AND Af <= 0.10)
		return 'core'
	# Surface criteria from paper: (5.7 < r < 7.0 AND Af >= 0.65) OR (r <= 5.7 AND Af >= 0.30)
	elif (5.7 < r < 7.0 and Af >= 0.65) or \
		 (r <= 5.7 and Af >= 0.30):
		return 'surface'
	# Boundary criteria from paper: (r >= 7.0 AND Af >= 0.30) OR (5.7 < r < 7.0 AND 0.10 < Af < 0.65) OR (r <= 5.7 AND Af <= 0.30)
	# These are the remaining conditions if not core or surface.
	elif (r >= 7.0 and Af >= 0.30) or \
		 (5.7 < r < 7.0 and 0.10 < Af < 0.65) or \
		 (r <= 5.7 and Af <= 0.30):
		return 'boundary'
	else:
		# This block should ideally not be reached if r and Af are valid numbers
		# and the above conditions correctly mirror the paper's mutually exclusive regions.
		# If 5.7 < r < 7.0 and Af is between 0.10 and 0.30 or between 0.30 and 0.65 this might need checking
		# Let's re-verify paper's conditions:
		# Core: (r>=7.0 & Af<0.30) or (5.7<r<7.0 & Af<=0.10)
		# Bnd: (r>=7.0 & Af>=0.30) or (5.7<r<7.0 & 0.10<Af<0.65) or (r<=5.7 & Af<=0.30)
		# Surf:(5.7<r<7.0 & Af>=0.65) or (r<=5.7 & Af>=0.30)
		
		# The order of checks matters. If it's not Core and not Surface, it should be Boundary by elimination if conditions are exhaustive.
		# The boundary condition (r >= 7.0 and Af >= 0.30) covers what's not core for r >= 7.0
		# The boundary condition (r <= 5.7 and Af <= 0.30) covers what's not surface for r <= 5.7
		# The boundary condition (5.7 < r < 7.0 and 0.10 < Af < 0.65) covers the middle ground for r in (5.7, 7.0)
		# It seems the conditions are exhaustive.
		print(f"Warning: Residue group with r={r:.2f}, Af={Af:.2f} did not fall into a defined Wisz region. Defaulting to 'boundary'. Check logic if this occurs frequently.")
		return 'boundary' # Defaulting to boundary as it's often an intermediate state.


def _get_averaged_param(param_prefix, region_i, region_j, wisz_params):
	"""
	Gets a Wisz parameter, averaging if regions are different.
	Region can be 'core', 'boundary', 'surface'.
	Example param_prefix: "QQ_alpha_" or "Delta_HbQP_G_"
	"""
	# Convert region names to suffices used in wisz_params keys (C, B, S)
	region_map = {'core': 'C', 'boundary': 'B', 'surface': 'S'}
	
	suffix_i = region_map.get(region_i)
	suffix_j = region_map.get(region_j)

	if not suffix_i or not suffix_j:
		# Fallback or error if regions are not recognized
		# Using surface parameters as a (potentially poor) fallback.
		print(f"Warning: Unrecognized region(s) '{region_i}' or '{region_j}' for parameter averaging. Defaulting to surface.")
		suffix_i = suffix_i or 'S'
		suffix_j = suffix_j or 'S'

	param_i_key = f"{param_prefix}{suffix_i}"
	param_j_key = f"{param_prefix}{suffix_j}"

	param_i = wisz_params.get(param_i_key)
	param_j = wisz_params.get(param_j_key)

	if param_i is None or param_j is None:
		raise ValueError(f"Wisz parameter missing for keys: {param_i_key} or {param_j_key}")

	if region_i == region_j:
		return param_i
	else:
		return (param_i + param_j) / 2.0

# --- Dielectric Constant Functions ---

def calculate_solv_epsilon_for_desolvation(region_i, wisz_params):
	"""
	Returns the region-dependent protein desolvation dielectric (solv_epsilon(p_i))
	as per discussion related to Eq. 12.
	Note: For surface region, desolvation penalty is considered negligible by Wisz.
		  The actual energy term (Eq 12) handles this with (1 - Af_i).
		  This function returns the epsilon to be used in the (1/epsilon_protein - 1/epsilon_water) part.
	"""
	region_map = {'core': 'C', 'boundary': 'B', 'surface': 'S'} # Surface added for completeness
	suffix_i = region_map.get(region_i)

	if region_i == 'core':
		return wisz_params['solv_epsilon_C']
	elif region_i == 'boundary':
		return wisz_params['solv_epsilon_B']
	elif region_i == 'surface':
		# In the surface region, groups are exposed to solvent and the desolvation penalty
		# is negligible; $\Delta G_{solv,i}(p_{i})$ is therefore set to 0 in this region.
		# This means the (1/solv_epsilon(pi) - 1/epsilon_w) term should effectively become zero
		# or handled such that Delta_G_solv is zero.
		# For now, returning a high dielectric (like water) implies 1/solv_epsilon is small.
		# The calling energy function (Eq 12) will handle the (1-Af) term.
		# If Af is high for surface, (1-Af) is small.
		# Let the energy function set Delta_G_solv to 0 for surface.
		# Here we provide a value that, if used in Eq 12, would lead to a small contribution
		# before the (1-Af_i) factor, or the caller can check region_i == 'surface'.
		return wisz_params.get('DEFAULT_WATER_DIELECTRIC', 80.0) # Effectively, use water dielectric.
	else:
		print(f"Warning: Unrecognized region '{region_i}' for solv_epsilon. Defaulting to boundary value.")
		return wisz_params['solv_epsilon_B']


def calculate_QP_epsilon(region_i, region_j, distance_ij, wisz_params):
	"""
	Calculates the distance- and region-dependent dielectric for charge-polar (QP) interactions.
	Implements Eq. 14 from Wisz & Hellinga (2003).

	Args:
		region_i (str): Region of the first group ('core', 'boundary', 'surface').
		region_j (str): Region of the second group.
		distance_ij (float): Distance between the groups.
		wisz_params (dict): Dictionary of Wisz & Hellinga parameters.

	Returns:
		float: The calculated dielectric constant.
	"""
	alpha_qp = _get_averaged_param("QP_alpha_", region_i, region_j, wisz_params)
	lambda_qp = _get_averaged_param("QP_lambda_", region_i, region_j, wisz_params)

	if distance_ij < 1e-6: # Avoid issues with zero distance if it occurs
		return alpha_qp # At r=0, exp term is 0, so epsilon = 1.0 + alpha_qp*(1-1) = 1.0. Or, return a high value to make energy small.
					   # Paper implies it starts at 1.0 and increases. (1.0 - e^0) = 0.
					   # So for r_ij = 0, epsilon = 1.0.
					   # However, this case should not occur for inter-atomic distances.
					   # Let's assume distance_ij is always > 0 for actual interactions.

	epsilon = 1.0 + alpha_qp * (1.0 - math.exp(-lambda_qp * distance_ij))
	return epsilon


def calculate_QQ_epsilon(region_i, region_j, distance_ij, wisz_params):
	"""
	Calculates the distance- and region-dependent dielectric for simple charge-charge (QQ)
	Coulombic interactions. Implements Eq. 16 from Wisz & Hellinga (2003).

	Args:
		region_i (str): Region of the first group ('core', 'boundary', 'surface').
		region_j (str): Region of the second group.
		distance_ij (float): Distance between the groups.
		wisz_params (dict): Dictionary of Wisz & Hellinga parameters.

	Returns:
		float: The calculated dielectric constant.
	"""
	alpha_qq = _get_averaged_param("QQ_alpha_", region_i, region_j, wisz_params)
	lambda_qq = _get_averaged_param("QQ_lambda_", region_i, region_j, wisz_params)

	if distance_ij < 1e-6: # Similar to QP_epsilon
		return 1.0 

	epsilon = 1.0 + alpha_qq * (1.0 - math.exp(-lambda_qq * distance_ij))
	return epsilon


def get_IP_epsilon(region_i, region_j, wisz_params):
	"""
	Gets the region-dependent dielectric for ion-pair (IP) interactions.
	This is position-dependent, not distance-dependent like QQ/QP epsilon.

	Args:
		region_i (str): Region of the first group ('core', 'boundary', 'surface').
		region_j (str): Region of the second group.
		wisz_params (dict): Dictionary of Wisz & Hellinga parameters.

	Returns:
		float: The dielectric constant for ion pairs.
	"""
	# The IP_epsilon parameters in Table II are single values per region (C, B, S).
	# We average them if the two interacting groups are in different regions.
	return _get_averaged_param("IP_epsilon_", region_i, region_j, wisz_params)


# Solution pKa values for ionizable groups (from Wisz paper text or Table I)
# These might be needed if pH-dependent charges are used. For fixed charges, less so.
SOLUTION_PKA = {
	"ASP": 4.0, "GLU": 4.4, "TYR": 9.6, "LYS": 10.4,
	"ARG": 12.0, "CYS": 8.3, "NTERM": 7.5, "CTERM": 3.8,
	"HIS_ND1": 6.5, # For N-delta protonated tautomer (HSD-like)
	"HIS_NE2": 6.9, # For N-epsilon protonated tautomer (HSE-like)
	# Heme propionate mentioned in paper (pKa ~4.4) - add if needed
}

def get_charge_at_ph(group_type, ph, base_charge_deprotonated, pka_solution):
	"""
	Calculates the average charge of an ionizable group at a given pH.
	Assumes single proton exchange.
	base_charge_deprotonated: -1 for acids, 0 for bases.
	"""
	# For acids (ASP, GLU, CYS, TYR, CTERM):
	# HA <-> H+ + A-
	# charge = base_charge_deprotonated + (1 - 1 / (1 + 10**(ph - pka_solution)))
	# Deprotonated fraction = 1 / (1 + 10**(pka - pH))
	# Protonated fraction = 1 / (1 + 10**(pH - pka))
	
	if base_charge_deprotonated < 0: # Acid
		# Charge = (deprotonated_fraction * -1) + (protonated_fraction * 0)
		charge = base_charge_deprotonated * (1 / (1 + 10**(pka_solution - ph)))
	else: # Base (LYS, ARG, NTERM, HIS)
		# B + H+ <-> BH+
		# Charge = (protonated_fraction * +1) + (deprotonated_fraction * 0)
		charge = 1.0 * (1 / (1 + 10**(ph - pka_solution)))
	return charge

def get_formal_charge_for_wisz(atom, residue_record, ph, solution_pka_map):
	"""
	Determines the formal charge of an ionizable group represented by 'atom' at a given pH.
	This is a more detailed placeholder.
	Returns the effective formal charge of the group (+1, -1, 0, or fractional).
	"""
	resn = residue_record.resn.upper()
	atom_name = atom.name.strip().upper()

	# Simplified examples - this needs to be robust and cover all Wisz ionizable types
	if resn == "ASP" and atom_name in ["OD1", "OD2", "CG"]: # CG included as part of group
		return get_charge_at_ph(resn, ph, -1.0, solution_pka_map.get(resn, 4.0))
	elif resn == "GLU" and atom_name in ["OE1", "OE2", "CD"]:
		return get_charge_at_ph(resn, ph, -1.0, solution_pka_map.get(resn, 4.4))
	elif resn == "LYS" and atom_name == "NZ":
		return get_charge_at_ph(resn, ph, 0.0, solution_pka_map.get(resn, 10.4)) # Base, deprotonated is 0
	elif resn == "ARG" and atom_name in ["NE", "CZ", "NH1", "NH2"]:
		# Arg is tricky as charge is delocalized. Paper says charge on NH1, NH2.
		# For simplicity, if any of these atoms, consider the group charge.
		return get_charge_at_ph(resn, ph, 0.0, solution_pka_map.get(resn, 12.0))
	elif resn == "HIS":
		# Average charge of HIS depends on which tautomer is more stable / pKa values
		# For simplicity, if pH < ~6, assume +1, else 0. This is a rough approximation.
		# A better way: pKa of HIP is around 6-7.
		pka_his = SOLUTION_PKA.get("HIS_NE2", 6.9) # Use one of the pKas as representative
		if ph < pka_his: return 1.0
		return 0.0
	elif atom_name == "N" and residue_record.atoms[0].id == atom.id: # N-terminus
		# Distinguish NTERM_PRO if your REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED does
		pka_nterm = SOLUTION_PKA.get("NTERM_PRO" if resn=="PRO" else "NTERM", 7.5)
		return get_charge_at_ph("NTERM", ph, 0.0, pka_nterm)
	elif atom_name == "OXT" or (atom_name == "O" and "OXT" not in [a.name for a in residue_record.atoms]): # C-terminus
		return get_charge_at_ph("CTERM", ph, -1.0, SOLUTION_PKA.get("CTERM", 3.8))
	elif resn == "CYS" and atom_name == "SG":
		return get_charge_at_ph(resn, ph, -1.0, solution_pka_map.get(resn, 8.3))
	elif resn == "TYR" and atom_name == "OH":
		return get_charge_at_ph(resn, ph, -1.0, solution_pka_map.get(resn, 9.6))
	
	return 0.0 # Default if not a clearly defined ionizable center for these terms

def calculate_wisz_electrostatic_energy(atom1, atom2, interaction_type, distance_ij, hori_instance):
	if hori_instance.dielectric_method != 'wisz' or not hori_instance.wisz_params:
		return 0.0

	wisz_params = hori_instance.wisz_params
	res_rec1 = hori_instance.residues.get((atom1.chain, atom1.resi))
	res_rec2 = hori_instance.residues.get((atom2.chain, atom2.resi))

	if not res_rec1 or not res_rec2: return 0.0

	group_type1_str = _get_ionizable_group_type_for_sasa(res_rec1.resn.upper(), atom1.name.strip().upper(), hori_instance.ref_sasa_ionizable)
	group_type2_str = _get_ionizable_group_type_for_sasa(res_rec2.resn.upper(), atom2.name.strip().upper(), hori_instance.ref_sasa_ionizable)

	# Region classification is per group (chain, resi, group_type)
	region1 = hori_instance.region_classification.get((atom1.chain, atom1.resi, group_type1_str), 'surface') if group_type1_str else 'surface'
	region2 = hori_instance.region_classification.get((atom2.chain, atom2.resi, group_type2_str), 'surface') if group_type2_str else 'surface'

	# Salt screening (standard Debye-Huckel for IP and general QQ/QP, not applied to specific DeltaG for HB)
	_I = hori_instance.ionic_strength
	# Use DEFAULT_WATER_DIELECTRIC from wisz_parameters
	_eps_w = wisz_params.get('DEFAULT_WATER_DIELECTRIC', hori_instance.wisz_params.get('DEFAULT_WATER_DIELECTRIC', 80.0))
	_T = hori_instance.temperature
	kappa_val = 0.0
	if _eps_w * _T > 1e-9:
		kappa_val = 50.29 * math.sqrt(abs(_I) / (_eps_w * _T)) # A^-1 [cite: 106]

	standard_salt_screening_factor = math.exp(-kappa_val * distance_ij) if distance_ij > 1e-6 else 1.0

	# Use COULOMB_CONSTANT from wisz_parameters for kcal/mol
	COULOMB_CONSTANT = hori_instance.wisz_params.get('COULOMB_CONSTANT', 332.0637) 
	energy = 0.0

	q_formal1 = get_formal_charge_for_wisz(atom1, res_rec1, hori_instance.ph, SOLUTION_PKA)
	q_formal2 = get_formal_charge_for_wisz(atom2, res_rec2, hori_instance.ph, SOLUTION_PKA)

	if interaction_type.startswith('hbond'):
		is_a1_formally_charged = abs(q_formal1) > 0.01 # Check if atom1 belongs to a formally charged group
		is_a2_formally_charged = abs(q_formal2) > 0.01 # Check if atom2 belongs to a formally charged group

		if is_a1_formally_charged and is_a2_formally_charged: # Charge-Charge H-bond (HbQQ) [cite: 126]
			energy = _get_averaged_param("Delta_HbQQ_G_", region1, region2, wisz_params)
		elif is_a1_formally_charged or is_a2_formally_charged: # Charge-Polar H-bond (HbQP) [cite: 124]
			energy = _get_averaged_param("Delta_HbQP_G_", region1, region2, wisz_params)
		else: # Polar-Polar H-bond (e.g., backbone-backbone, sidechain-backbone polar)
			  # Wisz model does not explicitly give a DeltaG term for these.
			  # The paper's Eq. 13 is for charge $q_i$ (formal) and polar $q_j$ (partial). [cite: 118]
			  # For polar-polar, an interpretation is needed. Using a QP-like calculation with partial charges:
			epsilon_qp = calculate_QP_epsilon(region1, region2, distance_ij, wisz_params)
			if epsilon_qp > 1e-6 and distance_ij > 1e-6:
				# Using partial charges (atom.charge from PDB2PQR) for both partners
				energy = (COULOMB_CONSTANT * atom1.charge * atom2.charge) / (epsilon_qp * distance_ij)
				# The paper's specific salt screening for non-HB QP (Eq 13) is e^(kappa*r*Af_avg).
				# For simplicity here, applying standard screening to this interpreted polar-polar term.
				energy *= standard_salt_screening_factor 
			else:
				energy = 0.0

	elif interaction_type == 'salt_bridge':
		# Wisz: oppositely charged ionizable groups, within 4A, NOT H-bonded [cite: 82, 127]
		if abs(q_formal1 * q_formal2) > 0.01 and (q_formal1 * q_formal2 < 0): # Oppositely formally charged
			epsilon_ip = get_IP_epsilon(region1, region2, wisz_params) # Position-dependent dielectric for IP
			if epsilon_ip > 1e-6 and distance_ij > 1e-6:
				energy = (COULOMB_CONSTANT * q_formal1 * q_formal2) / (epsilon_ip * distance_ij)
				# Apply standard salt screening to ion pair energy
				energy *= standard_salt_screening_factor
			else:
				energy = 0.0
		# If classified as salt_bridge by Hori's criteria (partial charges) but doesn't meet Wisz formal charge criteria,
		# it might be treated as a general QQ interaction below, or energy remains 0 from this block.

	return energy * 4.184 # conversion to kJ/Mol