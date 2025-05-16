from .classify_interactions import get_ring_data, is_aromatic
from .math_utils import distance, dist_3d, angle_between_vectors
import math

def compute_interaction_energy(a1, a2, dist, itype, residues, atoms, bonds, amber_nonbonded, user_params):
	if itype.startswith('hbond'):
		return hbond_energy(a1, a2, atoms, bonds)
	elif itype == 'salt_bridge':
		return salt_bridge_energy(a1, a2, dist)
	elif itype == 'vdw':
		return lennard_jones(a1, a2, dist, amber_nonbonded)
	elif itype == 'pi_pi':
		ring1 = get_ring_data(a1.chain, a1.resi, residues)
		ring2 = get_ring_data(a2.chain, a2.resi, residues)
		return pi_pi_energy(ring1, ring2, user_params)
	elif itype == 'cation_pi':
		return cation_pi_energy(a1, a2, residues)
	return 0.0

def hbond_energy(a1, a2, atoms, bonds):
	"""
	If we can do specialized formula => specialized_hbond_energy
	else => simple_hbond_energy
	"""
	if a1.atom_type.startswith('H'):
		h_atom = a1
		donor_atom = atoms[bonds[a1.id][0]]  # heavy neighbor
		acceptor_atom = a2
	else:
		h_atom = a2
		donor_atom = atoms[bonds[a2.id][0]]
		acceptor_atom = a1

	# must be N/O
	if not (donor_atom.atom_type.startswith('N') or donor_atom.atom_type.startswith('O')):
		return simple_hbond_energy(h_atom, acceptor_atom)
	if not (acceptor_atom.atom_type.startswith('N') or acceptor_atom.atom_type.startswith('O')):
		return simple_hbond_energy(h_atom, acceptor_atom)

	return specialized_hbond_energy(donor_atom, h_atom, acceptor_atom, atoms, bonds)

def simple_hbond_energy(h_atom, acceptor_atom):
	d = distance(h_atom, acceptor_atom)
	if d < 1e-12:
		return 0.0
	k = 332.06*4.184
	eps = 4.0
	return (k * h_atom.charge * acceptor_atom.charge) / (eps * d)

def specialized_hbond_energy(donor_atom, h_atom, acceptor_atom, atoms, bonds):
	"""
	Attempt the 4-distances approach for your custom formula.
	Otherwise fallback.
	"""
	rON = distance(donor_atom, acceptor_atom)
	rOH = distance(acceptor_atom, h_atom)

	c_accept = None
	for nbr_id in bonds[acceptor_atom.id]:
		nbr = atoms[nbr_id]
		if nbr.atom_type.startswith('C'):
			c_accept = nbr
			break
	if c_accept is None:
		return simple_hbond_energy(h_atom, acceptor_atom)
	rCH = distance(c_accept, h_atom)
	rCN = distance(c_accept, donor_atom)
	if min(rON, rOH, rCH, rCN) < 1e-6:
		return simple_hbond_energy(h_atom, acceptor_atom)

	factor = 332.06 * 4.184
	qd = donor_atom.charge
	qa = acceptor_atom.charge
	val = qd * qa * (1.0/rON + 1.0/rCH - 1.0/rOH - 1.0/rCN) * factor
	return val

def salt_bridge_energy(a1, a2, dist):
	k = 332.06 * 4.184  # kJ/mol
	eps = 4.0
	return (k * a1.charge * a2.charge) / (eps * dist)

def lennard_jones(a1, a2, dist, amber_nonbonded):
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

def cation_pi_energy(a1, a2, residues):
	"""
	Compute cation-π interaction energy.
	"""
	# Determine which atom is the cation and which is aromatic
	if a1.charge > 0.5 and is_aromatic(a2.resn):
		cation, aromatic = a1, a2
	elif a2.charge > 0.5 and is_aromatic(a1.resn):
		cation, aromatic = a2, a1
	else:
		return 0.0  # Not a valid cation-π interaction
	
	#get ring data 
	ring = get_ring_data(aromatic.chain, aromatic.resi, residues)
	if ring is None:
		return 0.0  # No valid ring data		
	# get centroid-cation distance
	cc_dist = dist_3d(ring['centroid'], (cation.x, cation.y, cation.z))
	
	# Parameters
	k_electrostatic = 332.0  # Coulomb constant in kJ/mol·Å·e
	optimal_distance = 4.0   # Å
	epsilon = 4.0
	# Electrostatic energy term
	energy = -1*(k_electrostatic * cation.charge) / (epsilon * cc_dist)

	# Dampen energy for non-optimal distances
	if cc_dist > optimal_distance:
		energy *= math.exp(-(cc_dist - optimal_distance) ** 2 / (2 * 1.5 ** 2))

	return energy