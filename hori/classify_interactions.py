from .math_utils import angle_between_vectors, cross_product, dist_3d
import math

def classify_interaction(distance_map, a1, a2, dist, atoms, bonds, residues, user_params):
	# 1) H-bond
	bool_is_h_bond, dist_h_bond = is_hbond(distance_map, a1, a2, atoms, bonds, user_params)
	if bool_is_h_bond:
		if dist_h_bond < 2.5:
			return 'hbond_strong'
		else:
			return 'hbond_weak'
	# 2) Salt bridge
	if dist <= user_params['salt_bridge_dist'] and (a1.charge * a2.charge < 0.0):
		return 'salt_bridge'
	# 3) pi-pi / cation-pi
	if is_aromatic(a1.resn) and is_aromatic(a2.resn):
		if is_pi_pi(a1, a2, residues, user_params):
			return 'pi_pi'
	elif is_cation_pi(a1, a2, residues, user_params):
		return 'cation_pi'
	# 4) vdw
	return 'vdw'

def is_hbond(distance_map, at1, at2, atoms, bonds, user_params):
	"""
	check if at1 and at2 form a hydrogen bond
	"""
	# which is hydrogen?
	if at1.atom_type.startswith('H'):
		h_atom, other = at1, at2
	elif at2.atom_type.startswith('H'):
		h_atom, other = at2, at1
	else:
		return False, None

	# acceptor must be O or N
	if not (other.atom_type.startswith('O') or other.atom_type.startswith('N')):
		return False, None

	# find donor heavy atom
	donors = bonds[h_atom.id]
	if len(donors) != 1:
		return False, None
	donor_atom = atoms[donors[0]]
	if donor_atom == other:
		return False, None #hydrogen cannot bond to back to the donor
	
	#donor atom check
	if not (donor_atom.atom_type.startswith('O') or donor_atom.atom_type.startswith('N')):
		return False, None #not h bond if donor is not O or N

	# distance check
	d_a_dist_key = (min(donor_atom.id, other.id), max(donor_atom.id, other.id))
	d_a_dist = distance_map.get(d_a_dist_key, None)
	h_a_dist_key = (min(h_atom.id, other.id), max(h_atom.id, other.id))
	h_a_dist = distance_map.get(h_a_dist_key, None)
	if d_a_dist is None or d_a_dist > user_params['d_a_dist'] or h_a_dist is None or h_a_dist > user_params['h_a_dist']:
		return False, None

	# angle check
	dh_vec = (h_atom.x - donor_atom.x, h_atom.y - donor_atom.y, h_atom.z - donor_atom.z)
	ah_vec = (h_atom.x - other.x, h_atom.y - other.y, h_atom.z - other.z)
	ang = angle_between_vectors(dh_vec, ah_vec)
	return (ang >= user_params['dha_angle']), d_a_dist

def is_aromatic(resn):
	return resn.upper() in ('PHE', 'TYR', 'TRP', 'HIS')

def is_pi_pi(a1, a2, residues, user_params):
	ring1 = get_ring_data(a1.chain, a1.resi, residues)
	ring2 = get_ring_data(a2.chain, a2.resi, residues)
	if (ring1 is None) or (ring2 is None):
		return False
	d = dist_3d(ring1['centroid'], ring2['centroid'])
	if d > user_params['pi_pi_dist']:
		return False
	angle = angle_between_vectors(ring1['normal'], ring2['normal'])
	# parallel or stacked
	if angle < user_params['pi_pi_angle'] or abs(180.0 - angle) < user_params['pi_pi_angle']:
		return True
	if (90.0-user_params['pi_pi_angle']) < angle < (90.0+user_params['pi_pi_angle']):
		return True 
	return False

def is_cation_pi(a1, a2, residues, user_params):
	# if dist > 6.0:
	# 	return False
	# check cation
	if a1.charge > 0.5 and is_aromatic(a2.resn):
		ring = get_ring_data(a2.chain, a2.resi, residues)
		if ring is None:
			return False
		d2 = dist_3d((a1.x, a1.y, a1.z), ring['centroid'])
		return (d2 < user_params['cation_pi_dist'])
	elif a2.charge > 0.5 and is_aromatic(a1.resn):
		ring = get_ring_data(a1.chain, a1.resi, residues)
		if ring is None:
			return False
		d2 = dist_3d((a2.x, a2.y, a2.z), ring['centroid'])
		return (d2 < user_params['cation_pi_dist'])
	return False

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