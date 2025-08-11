from .math_utils import angle_between_vectors, cross_product, dist_3d
import math
from .ramachandran import calculate_dihedral

def classify_interaction(distance_map, a1, a2, dist, atoms, bonds, residues, user_params):
	# 0) disulfide bond
	disulfide_details = is_disulfide_bond(a1, a2, dist, residues, user_params)
	if disulfide_details.pop('is_disulfide', False):
		disulfide_details['type'] = 'disulfide'
		return disulfide_details

	# 1) H-bond
	hbond_details = is_hbond(distance_map, a1, a2, atoms, bonds, user_params)
	if hbond_details.pop('is_hbond', False):
		if hbond_details['d_a_dist'] < 2.5:
			hbond_details['type'] = 'hbond_strong'
		else:
			hbond_details['type'] = 'hbond_weak'
		return hbond_details

	# 2) Salt bridge
	if dist <= user_params['salt_bridge_dist'] and (a1.charge * a2.charge < 0.0):
		return {'type': 'salt_bridge', 'distance': dist}

	# 3) pi-pi / cation-pi
	if is_aromatic(a1) and is_aromatic(a2):
		pi_pi_details = is_pi_pi(a1, a2, residues, user_params)
		if pi_pi_details.pop('is_pi_pi', False):
			pi_pi_details['type'] = 'pi_pi'
			return pi_pi_details
	else: # Use else to avoid re-checking aromaticity
		cation_pi_details = is_cation_pi(a1, a2, residues, user_params)
		if cation_pi_details.pop('is_cation_pi', False):
			cation_pi_details['type'] = 'cation_pi'
			return cation_pi_details
			
	# 4) vdw
	vdw_details = is_van_der_waals(a1, a2, dist, bonds, user_params)
	if vdw_details.pop('is_vdw', False):
		vdw_details['type'] = 'vdw'
		return vdw_details

	# none of the above
	return None

def is_van_der_waals(a1, a2, dist, bonds, user_params):
	"""
	Checks if the interaction between a1 and a2 qualifies as Van der Waals
	based on atomic radii and a tolerance. Excludes directly bonded atoms.

	Args:
		a1 (Atom): The first atom.
		a2 (Atom): The second atom.
		dist (float): The distance between a1 and a2.
		bonds (dict): A dictionary representing the bond graph {atom_id: [bonded_atom_ids,...]}.
		user_params (dict): Dictionary containing 'vdw_overlap_tolerance'.

	Returns:
		bool: True if it's a Van der Waals interaction, False otherwise.
	"""
	result = {'is_vdw': False}
	vdw_tol = user_params.get('vdw_overlap_tolerance', 0.5)

	# Ensure radii are valid numbers before using them
	if not (isinstance(a1.radius, (float, int)) and isinstance(a2.radius, (float, int))):
		return result # Cannot determine vdW interaction without valid radii

	sum_radii = a1.radius + a2.radius
	if dist < (sum_radii + vdw_tol):
		# Check that these atoms are not directly covalently bonded
		# (disulfide is handled separately; this is for other potential direct bonds)
		if bonds and a1.id in bonds and a2.id in bonds[a1.id]:
			return result # It's a direct covalent bond, not primarily a vdW interaction
		result['is_vdw'] = True
		if sum_radii > 0:
			result['rred'] = dist/sum_radii
		else:
			result['rred'] = 0
		return result
	return result

def _get_atom_from_record(residue_record, atom_name):
	"""Helper to get a specific atom from a ResidueRecord's atom list."""
	if residue_record:
		for atom in residue_record.atoms:
			if atom.name == atom_name:
				return atom
	return None

def _check_angle_within_options(angle, angle_options):
	"""Helper to check if an angle is within any of the (center, tolerance) options."""
	# Normalize angle to be within -180 to 180 for consistent checks, esp. for 180/-180 cases
	angle = (angle + 180) % 360 - 180

	for center, tolerance in angle_options:
		# Normalize center as well if it could be outside -180 to 180
		center_norm = (center + 180) % 360 - 180

		lower_bound = center_norm - tolerance
		upper_bound = center_norm + tolerance

		# Handle periodicity for centers near +/-180, e.g. center=180, tol=20 means (160, 200) or (160, -160)
		if center_norm == 180.0 or center_norm == -180.0: # Special handling for trans
			if abs(abs(angle) - 180.0) <= tolerance:
				return True
		elif lower_bound <= angle <= upper_bound:
			return True
	return False

def is_disulfide_bond(sg1_atom, sg2_atom, dist, all_residues_dict, user_params):
	"""
	Check for disulfide bond based on distance and dihedral angles, return geometric metrics if it is.
	sg1_atom, sg2_atom: The two sulfur atoms.
	all_residues_dict: The main dictionary mapping (chain, resi) to ResidueRecord.
	"""
	result = {'is_disulfide': False}
	# 1. Basic atom and residue type check
	if not (sg1_atom.resn == 'CYS' and sg1_atom.name == 'SG' and \
			sg2_atom.resn == 'CYS' and sg2_atom.name == 'SG'):
		return result

	# 2. S-S Distance check
	min_dist = user_params.get('disulfide_min_dist', 1.7)
	max_dist = user_params.get('disulfide_max_dist', 2.4)
	if not (min_dist <= dist <= max_dist):
		return result

	# Retrieve ResidueRecords for both Cysteines
	res1_rec = all_residues_dict.get((sg1_atom.chain, sg1_atom.resi))
	res2_rec = all_residues_dict.get((sg2_atom.chain, sg2_atom.resi))

	if not res1_rec or not res2_rec:
		return result # Should not happen if sg atoms are valid

	# Get necessary atoms for Cys1
	c1_cb = _get_atom_from_record(res1_rec, 'CB')
	c1_ca = _get_atom_from_record(res1_rec, 'CA')
	c1_n = _get_atom_from_record(res1_rec, 'N')

	# Get necessary atoms for Cys2
	c2_cb = _get_atom_from_record(res2_rec, 'CB')
	c2_ca = _get_atom_from_record(res2_rec, 'CA')
	c2_n = _get_atom_from_record(res2_rec, 'N')

	if not (c1_cb and c1_ca and c1_n and c2_cb and c2_ca and c2_n):
		return result # Cannot calculate dihedrals

	# Prepare coordinates for dihedral calculation
	p_c1_n  = (c1_n.x, c1_n.y, c1_n.z)
	p_c1_ca = (c1_ca.x, c1_ca.y, c1_ca.z)
	p_c1_cb = (c1_cb.x, c1_cb.y, c1_cb.z)
	p_sg1   = (sg1_atom.x, sg1_atom.y, sg1_atom.z)

	p_c2_n  = (c2_n.x, c2_n.y, c2_n.z)
	p_c2_ca = (c2_ca.x, c2_ca.y, c2_ca.z)
	p_c2_cb = (c2_cb.x, c2_cb.y, c2_cb.z)
	p_sg2   = (sg2_atom.x, sg2_atom.y, sg2_atom.z)

	# 3. chi-SS dihedral angle (CB1-SG1-SG2-CB2)
	chi_ss_val = calculate_dihedral(p_c1_cb, p_sg1, p_sg2, p_c2_cb)
	chi_ss_options = user_params.get('chi_ss_angle_opts', [(97.0, 30.0), (-87.0, 30.0)])
	if not _check_angle_within_options(chi_ss_val, chi_ss_options):
		return result

	# 4. chi1 dihedral angle (N-CA-CB-SG) for Cys1
	chi1_cys1_val = calculate_dihedral(p_c1_n, p_c1_ca, p_c1_cb, p_sg1)
	chi1_options = user_params.get('chi1_angle_opts', [(-60.0, 35.0), (60.0, 35.0), (180.0, 35.0)])
	if not _check_angle_within_options(chi1_cys1_val, chi1_options):
		return result

	# 5. chi1 dihedral angle (N-CA-CB-SG) for Cys2
	chi1_cys2_val = calculate_dihedral(p_c2_n, p_c2_ca, p_c2_cb, p_sg2)
	if not _check_angle_within_options(chi1_cys2_val, chi1_options):
		return result
	
	# If all checks passed, it's a disulfide bond
	result['is_disulfide'] = True
	result['dihedral_chi_ss'] = chi_ss_val
	result['dihedral_chi1_cys1'] = chi1_cys1_val
	result['dihedral_chi1_cys2'] = chi1_cys2_val

	return result

def is_hbond(distance_map, at1, at2, atoms, bonds, user_params):
	"""
	check if at1 and at2 form a hydrogen bond and returns metric data if they do.
	"""

	result = {'is_hbond': False}
	# which is hydrogen?
	if at1.atom_type.startswith('H'):
		h_atom, other = at1, at2
	elif at2.atom_type.startswith('H'):
		h_atom, other = at2, at1
	else:
		return result

	# acceptor must be O or N
	if not (other.atom_type.startswith('O') or other.atom_type.startswith('N')):
		return result

	# find donor heavy atom
	donors = bonds[h_atom.id]
	if len(donors) != 1:
		return result
	donor_atom = atoms[donors[0]]
	if donor_atom == other:
		return result #hydrogen cannot bond to back to the donor
	
	#donor atom check
	if not (donor_atom.atom_type.startswith('O') or donor_atom.atom_type.startswith('N')):
		return result #not h bond if donor is not O or N

	# distance check
	d_a_dist_key = (min(donor_atom.id, other.id), max(donor_atom.id, other.id))
	d_a_dist = distance_map.get(d_a_dist_key, None)
	h_a_dist_key = (min(h_atom.id, other.id), max(h_atom.id, other.id))
	h_a_dist = distance_map.get(h_a_dist_key, None)
	if d_a_dist is None or d_a_dist > user_params['d_a_dist'] or h_a_dist is None or h_a_dist > user_params['h_a_dist']:
		return result # not a hydrogen bond if distances are too large

	# angle check
	dh_vec = (h_atom.x - donor_atom.x, h_atom.y - donor_atom.y, h_atom.z - donor_atom.z)
	ah_vec = (h_atom.x - other.x, h_atom.y - other.y, h_atom.z - other.z)
	ang = angle_between_vectors(dh_vec, ah_vec)

	if ang >= user_params['dha_angle']:
		result['is_hbond'] = True
		result['d_a_dist'] = d_a_dist
		result['dha_angle'] = ang
		return result
	else:
		return result

def is_aromatic(atom):
	""" Check if the atom is one of the aromatic rings atom."""
	RING_DEFS = {'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
		'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
		'HIS': ['CG','ND1','CD2','CE1','NE2'],
		'TRP': ['CG','CD1','NE1','CE2','CD2','CE3','CZ2','CH2','CZ3']}
	resn = atom.resn.upper()
	if resn not in RING_DEFS:
		return False
	return atom.name.upper() in RING_DEFS[resn]

def is_pi_pi(a1, a2, residues, user_params):
	""" Check if a1 and a2 are aromatic residues forming a pi-pi interaction."""
	result = {'is_pi_pi': False}
	ring1 = get_ring_data(a1.chain, a1.resi, residues)
	ring2 = get_ring_data(a2.chain, a2.resi, residues)
	if (ring1 is None) or (ring2 is None):
		return result
	d = dist_3d(ring1['centroid'], ring2['centroid'])
	if d > user_params['pi_pi_dist']:
		return result
	angle = angle_between_vectors(ring1['normal'], ring2['normal'])
	# parallel or stacked
	is_parallel = angle < user_params['pi_pi_angle'] or abs(180.0 - angle) < user_params['pi_pi_angle']
	is_t_shaped = (90.0 - user_params['pi_pi_angle']) < angle < (90.0 + user_params['pi_pi_angle'])
	if is_parallel or is_t_shaped:
		result['is_pi_pi'] = True
		result['centroid_dist'] = d
		result['plane_angle'] = angle
		return result
	# not a pi-pi interaction
	return result

def is_cation_pi(a1, a2, residues, user_params):
	"""
	Check if a cation (positively charged atom) is interacting with an aromatic ring.
	"""
	result = {'is_cation_pi': False}

	cation_atom, aromatic_res_atom = (a1, a2) if a1.charge > 0.5 and is_aromatic(a2) else \
									 (a2, a1) if a2.charge > 0.5 and is_aromatic(a1) else \
									 (None, None)
	if not cation_atom:
		return result
	
	ring = get_ring_data(aromatic_res_atom.chain, aromatic_res_atom.resi, residues)
	if ring is None:
		return result
	
	cation_coords = (cation_atom.x, cation_atom.y, cation_atom.z)
	dist_to_centroid = dist_3d(cation_coords, ring['centroid'])

	if dist_to_centroid < user_params['cation_pi_dist']:
		result['is_cation_pi'] = True
		result['cation_centroid_dist'] = dist_to_centroid
		return result
	
	return result

def get_ring_data(chain, resi, residues):
	"""
	Return centroid & normal for ring in PHE/TYR/HIS/TRP if found.
	"""
	ring_defs = {
		'PHE': ['CG','CD1','CD2','CE1','CE2','CZ'],
		'TYR': ['CG','CD1','CD2','CE1','CE2','CZ'],
		'HIS': ['CG','ND1','CD2','CE1','NE2'],
		'TRP': ['CG','CD1','NE1','CE2','CD2','CE3','CZ2','CH2','CZ3']
	}
	rec = residues.get((chain, resi))
	if rec is None:
		return None
	resn = rec.resn.upper()
	if resn not in ring_defs:
		return None
	ring_names = ring_defs[resn]
	coords = []
	for a in rec.atoms:
		if a.name in ring_names:
			coords.append((a.x, a.y, a.z))
	if len(coords) < 3:
		return None

	cx = sum(c[0] for c in coords) / len(coords)
	cy = sum(c[1] for c in coords) / len(coords)
	cz = sum(c[2] for c in coords) / len(coords)
	centroid = (cx, cy, cz)

	v1 = (coords[0][0] - cx, coords[0][1] - cy, coords[0][2] - cz)
	v2 = (coords[1][0] - cx, coords[1][1] - cy, coords[1][2] - cz)
	normal = cross_product(v1, v2)
	mag = math.sqrt(normal[0]**2 + normal[1]**2 + normal[2]**2)
	if mag < 1e-6:
		return None
	normal = (normal[0]/mag, normal[1]/mag, normal[2]/mag)
	return {'centroid': centroid, 'normal': normal}