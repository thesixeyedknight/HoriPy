import math
from .spatial_utils import get_atoms_in_cutoff_from_grid

def _calculate_local_polarizability_density(target_atom, all_atoms_dict, grid_data, pike_params):
	"""
	Calculates the local polarizability density (rho_alpha) around a target atom.
	rho_alpha = (Sum of polarizabilities of all protein atoms and waters in the sphere) / (Volume of the sphere)

	Args:
		target_atom (Atom): The Atom object for which rho_alpha is calculated.
		all_atoms_dict (dict): Main Hori dictionary mapping atom ID to Atom object.
		grid_data (dict): Spatial grid data.
		pike_params (dict): Pike & Nanda parameters, including:
							- 'AMBER_TO_PIKE_ATOM_TYPE_MAP'
							- 'PIKE_ATOMIC_POLARIZABILITIES'
							- 'PIKE_ATOMIC_VOLUMES'
							- 'POLARIZABILITY_INFLUENCE_RADIUS' (should be 9.0 Å)

	Returns:
		float: The local polarizability density (dimensionless).
	"""
	if not pike_params:
		print("Warning: Pike & Nanda parameters missing in _calculate_local_polarizability_density.")
		return 0.0

	amber_to_pike_map = pike_params.get('AMBER_TO_PIKE_ATOM_TYPE_MAP', {})
	atomic_polarizabilities = pike_params.get('PIKE_ATOMIC_POLARIZABILITIES', {})
	atomic_volumes = pike_params.get('PIKE_ATOMIC_VOLUMES', {})
	influence_radius = pike_params.get('POLARIZABILITY_INFLUENCE_RADIUS')
	alpha_HOH = atomic_polarizabilities.get('HOH') # Polarizability of one water molecule
	V_HOH = atomic_volumes.get('HOH')             # Volume of one water molecule

	if influence_radius is None or influence_radius <= 0:
		print("Warning: Invalid POLARIZABILITY_INFLUENCE_RADIUS in Pike & Nanda parameters.")
		return 0.0
	if not amber_to_pike_map:
		print("Warning: AMBER_TO_PIKE_ATOM_TYPE_MAP missing in Pike & Nanda parameters.")
		return 0.0

	# 1. Get protein neighbors within the influence radius
	protein_neighbors = get_atoms_in_cutoff_from_grid(
		(target_atom.x, target_atom.y, target_atom.z),
		target_atom.id,
		all_atoms_dict,
		grid_data,
		influence_radius
	)

	# 2. Calculate sum of polarizabilities and volume occupied by these protein neighbors
	sum_alpha_protein_in_sphere = 0.0
	volume_occupied_by_protein_in_sphere = 0.0

	for neighbor_atom in protein_neighbors:
		amber_type = neighbor_atom.atom_type
		pike_type = amber_to_pike_map.get(amber_type)

		if pike_type:
			alpha_p = atomic_polarizabilities.get(pike_type)
			V_p = atomic_volumes.get(pike_type)

			if alpha_p is not None:
				sum_alpha_protein_in_sphere += alpha_p
			else:
				print(f"Debug: Missing polarizability for Pike type {pike_type} (AMBER: {amber_type})")


			if V_p is not None:
				volume_occupied_by_protein_in_sphere += V_p
			else:
				print(f"Debug: Missing volume for Pike type {pike_type} (AMBER: {amber_type})")


	# 3. Calculate volume of the influence sphere
	V_sphere = (4.0 / 3.0) * math.pi * (influence_radius ** 3)
	if V_sphere < 1e-9: # Avoid division by zero if radius is tiny
		return 0.0

	# 4. Calculate contribution from water
	sum_alpha_water_in_sphere = 0.0
	if alpha_HOH is not None and V_HOH is not None and V_HOH > 1e-9:
		volume_for_water = max(0.0, V_sphere - volume_occupied_by_protein_in_sphere)
		num_effective_waters = volume_for_water / V_HOH
		sum_alpha_water_in_sphere = num_effective_waters * alpha_HOH
	# else:
		# print("Debug: Water parameters (alpha_HOH or V_HOH) missing or invalid. No water contribution.")


	# 5. Calculate total sum of polarizabilities in the sphere
	Total_Sum_Alpha_in_Sphere = sum_alpha_protein_in_sphere + sum_alpha_water_in_sphere

	# 6. Calculate polarizability density (rho_alpha)
	polarizability_density = Total_Sum_Alpha_in_Sphere / V_sphere
	
	return polarizability_density

def _convert_polarizability_density_to_epsilon(polarizability_density, pike_params):
	"""
	Converts polarizability density (rho_alpha) to local dielectric constant (epsilon_local)
	using the formula: epsilon = intercept + scale_factor * rho_alpha.

	Args:
		polarizability_density (float): The dimensionless local polarizability density.
		pike_params (dict): Pike & Nanda parameters, including 'EPSILON_FUNCTION_PARAMS'.

	Returns:
		float: The calculated local dielectric constant.
	"""
	if not pike_params:
		print("Warning: Pike & Nanda parameters missing in _convert_polarizability_density_to_epsilon.")
		return 1.0

	epsilon_func_params = pike_params.get('EPSILON_FUNCTION_PARAMS', {})
	intercept = epsilon_func_params.get('intercept', 1.0)
	scale_factor = epsilon_func_params.get('scale_factor_for_N_alpha_product', 4.0 * math.pi)

	epsilon_local = intercept + scale_factor * polarizability_density
	
	return max(1.0, epsilon_local) # Ensure dielectric is physically realistic (>=1)

def calculate_pike_effective_dielectric_for_interaction(atom1, atom2, hori_instance, pike_params, grid_data):
	"""
	Calculates the effective dielectric constant for an interaction between atom1 and atom2
	using the Pike & Nanda model. It averages the local dielectric constants at the
	positions of atom1 and atom2.

	Args:
		atom1 (Atom): The first Atom object.
		atom2 (Atom): The second Atom object.
		hori_instance: The main Hori instance (to access hori_instance.atoms).
		pike_params (dict): Pike & Nanda parameters.
		grid_data (dict): Spatial grid data.

	Returns:
		float: The effective dielectric constant for the interaction.
	"""
	# Default to a high dielectric (like bulk water) if prerequisites are missing
	default_epsilon = pike_params.get('BULK_SOLVENT_DIELECTRIC_EPSILON_S', 80.0)

	if not hori_instance or not pike_params or not grid_data:
		print("Warning: Missing inputs for Pike & Nanda effective dielectric calculation. Returning default bulk.")
		return default_epsilon

	# Calculate local polarizability density and then epsilon for atom1
	rho_alpha1 = _calculate_local_polarizability_density(atom1, hori_instance.atoms, grid_data, pike_params)
	epsilon_local1 = _convert_polarizability_density_to_epsilon(rho_alpha1, pike_params)

	# Calculate local polarizability density and then epsilon for atom2
	rho_alpha2 = _calculate_local_polarizability_density(atom2, hori_instance.atoms, grid_data, pike_params)
	epsilon_local2 = _convert_polarizability_density_to_epsilon(rho_alpha2, pike_params)

	# Average the two local dielectric constants for the interaction path
	epsilon_effective = (epsilon_local1 + epsilon_local2) / 2.0
	
	return max(1.0, epsilon_effective) # Ensure final effective dielectric is at least 1.0