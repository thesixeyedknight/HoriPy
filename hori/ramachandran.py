import math

def calculate_dihedral(p0, p1, p2, p3):
	"""
	Calculate dihedral angle between four points.

	Args:
		p0, p1, p2, p3: Tuples representing the (x, y, z) coordinates of the four atoms.

	Returns:
		Dihedral angle in degrees.
	"""
	def vector(a, b):
		return (b[0] - a[0], b[1] - a[1], b[2] - a[2])

	def normalize(v):
		mag = math.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
		return (v[0]/mag, v[1]/mag, v[2]/mag) if mag > 0 else (0.0, 0.0, 0.0)

	def cross(v1, v2):
		return (
			v1[1]*v2[2] - v1[2]*v2[1],
			v1[2]*v2[0] - v1[0]*v2[2],
			v1[0]*v2[1] - v1[1]*v2[0]
		)

	def dot(v1, v2):
		return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]

	b0 = vector(p0, p1)
	b1 = vector(p1, p2)
	b2 = vector(p2, p3)

	# Normalize b1 so that it does not influence magnitude of vector rejections
	b1_norm = normalize(b1)

	# Compute cross products
	v = cross(b0, b1)
	w = cross(b1, b2)

	x = dot(v, w)
	y = dot(cross(v, b1_norm), w)

	angle = math.degrees(math.atan2(y, x))
	return angle

def get_atom(residue, atom_name):
	"""
	Retrieve an atom from a residue by its name.

	Args:
		residue: ResidueRecord object.
		atom_name: String name of the atom.

	Returns:
		Atom object if found, else None.
	"""
	for atom in residue.atoms:
		if atom.name == atom_name:
			return atom
	return None

def compute_ramachandran_angles(residues):
	"""
	Compute the Ramachandran angles (phi, psi, omega) for each residue.

	Args:
		residues: Dictionary with keys as (chain, resi) and values as ResidueRecord.

	Returns:
		Dictionary with the same keys as input, where each ResidueRecord has phi, psi, and omega populated.
	"""
	# Sort residues by chain and residue index
	sorted_keys = sorted(residues.keys(), key=lambda x: (x[0], x[1]))
	sorted_residues = [residues[key] for key in sorted_keys]

	for i, res in enumerate(sorted_residues):
		# Initialize angles as None
		phi = None
		psi = None
		omega = None

		# Check if backbone atoms are present in current residue
		N = get_atom(res, 'N')
		CA = get_atom(res, 'CA')
		C = get_atom(res, 'C')

		if not (N and CA and C):
			# Cannot compute angles without backbone atoms
			sorted_residues[i] = res._replace(phi=phi, psi=psi, omega=omega)
			continue

		# Compute Phi
		if i > 0:
			prev_res = sorted_residues[i - 1]
			C_prev = get_atom(prev_res, 'C')
			if C_prev:
				phi = calculate_dihedral(
					(C_prev.x, C_prev.y, C_prev.z),
					(N.x, N.y, N.z),
					(CA.x, CA.y, CA.z),
					(C.x, C.y, C.z)
				)

		# Compute Psi
		if i < len(sorted_residues) - 1:
			next_res = sorted_residues[i + 1]
			N_next = get_atom(next_res, 'N')
			if N_next:
				psi = calculate_dihedral(
					(N.x, N.y, N.z),
					(CA.x, CA.y, CA.z),
					(C.x, C.y, C.z),
					(N_next.x, N_next.y, N_next.z)
				)

		# Compute Omega
		if i < len(sorted_residues) - 1:
			next_res = sorted_residues[i + 1]
			N_next = get_atom(next_res, 'N')
			CA_next = get_atom(next_res, 'CA')
			if N_next and CA_next:
				omega = calculate_dihedral(
					(CA.x, CA.y, CA.z),
					(C.x, C.y, C.z),
					(N_next.x, N_next.y, N_next.z),
					(CA_next.x, CA_next.y, CA_next.z)
				)

		# Update the residue record
		sorted_residues[i] = res._replace(phi=phi, psi=psi, omega=omega)

	# Reconstruct the residues dictionary
	updated_residues = {key: res for key, res in zip(sorted_keys, sorted_residues)}
	return updated_residues
