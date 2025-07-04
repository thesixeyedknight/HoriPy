from .classify_interactions import get_ring_data, is_aromatic
from .math_utils import distance, dist_3d, angle_between_vectors
from .pike_nanda_utils import calculate_pike_effective_dielectric_for_interaction 
import math


def compute_interaction_energy(a1, a2, dist, itype, residues, atoms, bonds, amber_nonbonded, user_params, hori_instance=None):
	"""
	hori_instance is required if user_params['dielectric_method'] is 'pike_nanda'.
	"""
	# For 'bulk' or 'pike_nanda', the specific energy functions below will handle dielectric.
	# Pass user_params and hori_instance to them.

	if itype.startswith('hbond'):
		# Pass 'dist' which is likely the D-A or H-A distance map value.
		# hbond_energy will need to recalculate specific distances like H-A if 'dist' is D-A.
		# Or, ensure 'dist' passed here is relevant (e.g., H-A if that's what simple_hbond_energy uses primarily)
		# For now, we pass `dist` as the general interaction distance from distance_map.
		return hbond_energy(a1, a2, dist, atoms, bonds, user_params, hori_instance)
	elif itype == 'salt_bridge':
		return salt_bridge_energy(a1, a2, dist, user_params, hori_instance)
	elif itype == 'vdw':
		return lennard_jones(a1, a2, dist, amber_nonbonded) # Unaffected by these dielectric models
	elif itype == 'pi_pi':
		ring1 = get_ring_data(a1.chain, a1.resi, residues)
		ring2 = get_ring_data(a2.chain, a2.resi, residues)
		return pi_pi_energy(ring1, ring2, user_params) # Generally not using Coulombic dielectric directly
	elif itype == 'cation_pi':
		# 'dist' here is the atom-atom distance from distance_map.
		# cation_pi_energy calculates centroid-cation distance internally.
		return cation_pi_energy(a1, a2, dist, residues, user_params, hori_instance)
	elif itype == 'disulfide':
		return -251.0  # Disulfide bond energy is constant in this context
	return 0.0

# Helper: _get_hbond_dielectric (NEW)
def _get_hbond_dielectric(donor_heavy_atom, acceptor_heavy_atom, user_params, hori_instance):
	epsilon = user_params.get('bulk_dielectric_value', 4.0) # Start with bulk or default

	if user_params.get('dielectric_method') == 'pike_nanda':
		pike_params = user_params.get('pike_nanda_params')
		spatial_grid_data = user_params.get('spatial_grid_data')
		if hori_instance and pike_params and spatial_grid_data:
			eff_eps = calculate_pike_effective_dielectric_for_interaction(
				donor_heavy_atom, acceptor_heavy_atom, hori_instance, pike_params, spatial_grid_data
			)
			epsilon = eff_eps
		else:
			# This warning might be redundant if analysis.py already fell back user_params['dielectric_method']
			# but good for direct debugging here.
			print(f"Warning: Pike & Nanda prerequisites missing for H-bond dielectric. Using epsilon={epsilon}.")
	# 'bulk' is handled by the initial value of epsilon from user_params.
	# 'wisz' is handled separately in compute_interaction_energy.
	return max(1.0, epsilon)

def hbond_energy(a1, a2, dist_map_value, atoms, bonds, user_params, hori_instance=None): # Added dist_map_value
	# 1. Identify Donor (D), Hydrogen (H), Acceptor (A)
	if a1.atom_type.startswith('H'):
		h_atom = a1
		acceptor_heavy_atom = a2
	elif a2.atom_type.startswith('H'):
		h_atom = a2
		acceptor_heavy_atom = a1
	else:
		# This case should ideally not be reached if 'itype' was correctly 'hbond'
		return 0.0 

	donor_heavy_atom_ids = bonds.get(h_atom.id, [])
	if not donor_heavy_atom_ids: return 0.0 # Hydrogen not bonded to anything
	donor_heavy_atom = atoms[donor_heavy_atom_ids[0]]

	# Check if donor and acceptor are N/O for specialized, otherwise use simple
	is_special_type = (donor_heavy_atom.atom_type.startswith('N') or donor_heavy_atom.atom_type.startswith('O')) and \
					  (acceptor_heavy_atom.atom_type.startswith('N') or acceptor_heavy_atom.atom_type.startswith('O'))

	# 2. Determine Epsilon (common for both simple and specialized formulas)
	# The dielectric is primarily influenced by the environment of the heavy atoms.
	epsilon = _get_hbond_dielectric(donor_heavy_atom, acceptor_heavy_atom, user_params, hori_instance)
	
	# 3. Calculate energy using the appropriate formula
	if is_special_type:
		# Find C neighbor of acceptor for specialized formula
		c_acceptor_neighbor = None
		if acceptor_heavy_atom.id in bonds:
			for nbr_id in bonds[acceptor_heavy_atom.id]:
				nbr = atoms.get(nbr_id)
				if nbr and nbr.atom_type.startswith('C'): # Assuming AMBER 'C', 'CT', 'CA' etc.
					c_acceptor_neighbor = nbr
					break
		return _specialized_hbond_formula(donor_heavy_atom, h_atom, acceptor_heavy_atom, c_acceptor_neighbor, epsilon)
	else:
		dist_ha = distance(h_atom, acceptor_heavy_atom) # Calculate H-A distance specifically
		return _simple_hbond_formula(h_atom, acceptor_heavy_atom, dist_ha, epsilon)


def _simple_hbond_formula(h_atom, acceptor_atom, dist_ha, epsilon): # Now takes epsilon
	# dist_ha is the H-A distance
	if dist_ha < 1e-12:
		return 0.0
	k_coulomb = 332.06 * 4.184 # kJ/mol A e^-2
	# Energy calculation using partial charges of H and Acceptor
	return (k_coulomb * h_atom.charge * acceptor_atom.charge) / (epsilon * dist_ha)


def _specialized_hbond_formula(donor_atom, h_atom, acceptor_atom, c_acceptor_neighbor, epsilon): # Now takes epsilon
	rON = distance(donor_atom, acceptor_atom)
	rOH = distance(acceptor_atom, h_atom) # H-A distance
	
	if c_acceptor_neighbor is None: # Fallback if C neighbor of acceptor not found
		# print("Debug: C_acceptor_neighbor not found for specialized H-bond, falling back to simple formula logic with H-A.")
		return _simple_hbond_formula(h_atom, acceptor_atom, rOH, epsilon)

	rCH = distance(c_acceptor_neighbor, h_atom)
	rCN = distance(c_acceptor_neighbor, donor_atom)

	if min(rON, rOH, rCH, rCN) < 1e-6: # Avoid division by zero if any critical distance is ~0
		# Fallback to simple formula if distances are problematic
		# print("Debug: Problematic distance in specialized H-bond, falling back to simple formula logic with H-A.")
		return _simple_hbond_formula(h_atom, acceptor_atom, rOH, epsilon)

	# The original formula: qd * qa * (1.0/rON + 1.0/rCH - 1.0/rOH - 1.0/rCN) * factor
	# Here qd is donor_atom.charge and qa is acceptor_atom.charge. Factor includes epsilon.
	factor_with_eps = (332.06 * 4.184) / epsilon
	val = donor_atom.charge * acceptor_atom.charge * \
		  (1.0/rON + 1.0/rCH - 1.0/rOH - 1.0/rCN) * factor_with_eps
	return val

def salt_bridge_energy(a1, a2, dist, user_params, hori_instance=None): # Added user_params, hori_instance
	k_coulomb = 332.06 * 4.184  # kJ/mol A e^-2
	epsilon = user_params.get('bulk_dielectric_value', 4.0) # Default to bulk or a fallback
	if user_params.get('dielectric_method') == 'pike_nanda':
		pike_params = user_params.get('pike_nanda_params')
		spatial_grid_data = user_params.get('spatial_grid_data')
		if hori_instance and pike_params and spatial_grid_data:
			# For salt bridge, a1 and a2 are the charged atoms defining the interaction
			eff_eps = calculate_pike_effective_dielectric_for_interaction(
				a1, a2, hori_instance, pike_params, spatial_grid_data
			)
			epsilon = eff_eps
		else:
			# analysis.py should have set user_params['dielectric_method'] to 'bulk' if prerequisites failed.
			# This warning is an additional safeguard or for debugging.
			print(f"Warning: Pike & Nanda prerequisites missing for salt_bridge ({a1.id}-{a2.id}). Using epsilon={epsilon}.")
			# Epsilon remains the bulk value if Pike & Nanda setup failed.

	# Note: 'bulk' case is covered by the initial assignment of epsilon using user_params.get('bulk_dielectric_value',...)
	if dist < 1e-6: return 0.0 # Avoid division by zero
	return (k_coulomb * a1.charge * a2.charge) / (max(1.0, epsilon) * dist)

def lennard_jones(a1, a2, dist, amber_nonbonded):
	if dist < 1e-6:
		return 0.0
	s1, e1 = amber_nonbonded.get(a1.atom_type, (0.34, 0.2))
	s2, e2 = amber_nonbonded.get(a2.atom_type, (0.34, 0.2))
	sigma = 0.5*(s1 + s2)
	epsilon = math.sqrt(e1*e2)
	sr = sigma/dist
	sr6 = sr**6
	sr12 = sr6**2
	return 4.0*epsilon*(sr12 - sr6)

def pi_pi_energy(ring1, ring2, user_params):
	"""
	Compute π-π interaction energy, distinguishing between parallel and T-stacking.
	"""
	d = dist_3d(ring1['centroid'], ring2['centroid'])
	angle = angle_between_vectors(ring1['normal'], ring2['normal'])

	# Parameters
	parallel_strength = -5.0  # kJ/mol
	t_stack_strength = -2.0   # kJ/mol
	optimal_distance = 3.5    # Å
	distance_tolerance = 0.75 # Å

	# Distance-dependent term
	d_term = math.exp(-((d - optimal_distance) ** 2) / (2 * distance_tolerance ** 2))

	# Parallel stacking
	if angle < user_params['pi_pi_angle'] or abs(180.0 - angle) < user_params['pi_pi_angle']:
		return parallel_strength * d_term

	# T-stacking
	if (90.0-user_params['pi_pi_angle']) < angle < (90.0+user_params['pi_pi_angle']):
		return t_stack_strength * d_term

	return 0.0  # No interaction

def cation_pi_energy(a1, a2, dist_atom_atom, residues, user_params, hori_instance=None): # Added dist_atom_atom, user_params, hori_instance
	# Determine which atom is the cation and which is representative of the aromatic ring
	if a1.charge > 0.5 and is_aromatic(a2.resn): # a1 is cation, a2 is on aromatic ring
		cation_atom, aromatic_ring_atom_ref = a1, a2
	elif a2.charge > 0.5 and is_aromatic(a1.resn): # a2 is cation, a1 is on aromatic ring
		cation_atom, aromatic_ring_atom_ref = a2, a1
	else:
		return 0.0 # Not a valid cation-pi pair based on initial check

	ring_data = get_ring_data(aromatic_ring_atom_ref.chain, aromatic_ring_atom_ref.resi, residues)
	if ring_data is None:
		return 0.0

	# Distance for energy formula is cation to ring centroid
	cc_dist = dist_3d(ring_data['centroid'], (cation_atom.x, cation_atom.y, cation_atom.z))
	if cc_dist < 1e-6: return 0.0

	epsilon = user_params.get('bulk_dielectric_value', 4.0) # Default to bulk or fallback

	if user_params.get('dielectric_method') == 'pike_nanda':
		pike_params = user_params.get('pike_nanda_params')
		spatial_grid_data = user_params.get('spatial_grid_data')
		if hori_instance and pike_params and spatial_grid_data:
			# For cation-pi, dielectric is between cation and the pi-system.
			# We use the initially identified aromatic_ring_atom_ref as the representative
			# atom from the ring for calculating its local dielectric environment.
			eff_eps = calculate_pike_effective_dielectric_for_interaction(
				cation_atom, aromatic_ring_atom_ref, hori_instance, pike_params, spatial_grid_data
			)
			epsilon = eff_eps
		else:
			print(f"Warning: Pike & Nanda prerequisites missing for cation-pi. Using epsilon={epsilon}.")
	
	k_electrostatic = 332.06 * 4.184 # kJ/mol A e^-2
	optimal_distance_cp = 4.0 # Å, for the damping term, not for dielectric calc

	# Standard simplified electrostatic model for cation-pi
	# Assumes an effective negative charge on the pi ring interacting with the cation.
	# The charge of the pi system is effectively -1 (if cation is +ve) or +1 (if cation is -ve) for attraction.
	effective_pi_charge_sign = -math.copysign(1.0, cation_atom.charge) 
	energy = (k_electrostatic * cation_atom.charge * effective_pi_charge_sign) / (max(1.0, epsilon) * cc_dist)

	# Dampen energy for non-optimal distances (empirical factor)
	if cc_dist > optimal_distance_cp:
		energy *= math.exp(-((cc_dist - optimal_distance_cp)**2) / (2 * 1.5**2)) # Damping factor

	return energy
