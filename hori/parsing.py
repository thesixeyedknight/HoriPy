# parsing.py
def test_using_number_of_columns(lines):
	"""
	Determine the file format (PDB or CIF) based on the number of columns in the ATOM lines.
	"""
	found_atom = False
	for line in lines:
		if line.startswith('ATOM'):
			cols = line.split()
			if len(cols) == 12:
				return 'pdb'
			elif len(cols) == 18 or len(cols) == 21:
				return 'cif'
			else:
				raise ValueError(f"Unknown file format: Unexpected column count {len(cols)} in line:\n{line}")
	if not found_atom:
		raise ValueError("No 'ATOM' lines found in the file; cannot determine file format")


def read_pdb(lines):
	"""
	Parse PDB file lines to extract ATOM records and chains.
	"""
	pdb_for_haad = []
	chains = set()
	for line in lines:
		if line.startswith('ATOM'):
			pdb_for_haad.append(line)
			chains.add(line[21])
		elif line.strip() == 'ENDMDL':
			break
	return pdb_for_haad, chains

def make_pdb_from_cif(cif_lines):
	pdb_lines = []
	for line in cif_lines:
		if line.startswith('ATOM'):
			cols = line.split()
			if cols[-1]=='1':
				pdb_line = (
					"ATOM  {:>5s} {:^4s} {:>3s} {:>1s}{:>4s}    "
					"{:>8s}{:>8s}{:>8s}{:>6s}{:>6s}          {:>2s}\n"
					.format(
						cols[1],
						cols[3].strip('"'), cols[5], cols[6], cols[8],
						cols[10], cols[11], cols[12],
						cols[13], cols[14], cols[2]
					)
				)
				pdb_lines.append(pdb_line)
	return pdb_lines

def read_cif(lines):
	"""
	Parse CIF file lines to extract ATOM records and chains.
	"""
	pdb_for_haad = []
	chains = set()
	for line in lines:
		if line.startswith('ATOM'):
			cols = line.split()
			if cols[-1] == '1':
				pdb_for_haad.append(line)
				chains.add(cols[6])
	return pdb_for_haad, chains


def read_file(filename):
	"""
	Read a structure file (PDB or CIF) and determine its type.
	"""
	with open(filename, 'r') as f:
		lines = f.readlines()

	if lines and lines[0].startswith('HEADER'):
		for line in lines:
			if line.startswith('ATOM'):
				file_type = 'pdb'
				break
	elif lines and lines[0].startswith('data_'):
		for line in lines:
			if line.startswith('ATOM'):
				file_type = 'cif'
				break
	else:
		file_type = test_using_number_of_columns(lines)

	if file_type == 'pdb':
		pdb_for_haad, chains = read_pdb(lines)
	elif file_type == 'cif':
		pdb_for_haad, chains = read_cif(lines)
	else:
		raise ValueError('Unknown file format')

	return file_type, pdb_for_haad, chains

pqr_to_amber = {
	# Common backbone heavy atoms
	'N':   'N',    # backbone N often is type 'N' or 'N3'
	'CA':  'CT',   # backbone alpha-carbon
	'C':   'C',    # backbone carbonyl carbon
	'O':   'O',    # backbone carbonyl oxygen
	'OXT': 'O2',   # terminal carboxyl oxygen (often -1 charge)
	'H':   'H1',   # backbone H on the amide nitrogen

	# Lys/Arg side-chain nitrogens often called 'N3' in AMBER
	'NZ':  'N3',   # Lys terminal N
	'NH1': 'N3',   # Arg
	'NH2': 'N3',   # Arg
	'NE':  'N3',   # Arg or Lys NE

	# Asp/Glu side-chain oxygens often 'O2'
	'OD1': 'O2',
	'OD2': 'O2',
	'OE1': 'O2',
	'OE2': 'O2',

	# Ser/Thr hydroxyl oxygen -> 'OH' or 'OS'
	'OG':  'OH',   # Ser
	'OG1': 'OH',   # Thr

	# Cys
	'SG':  'S',    
	# Met
	'SD':  'S',    

	# Some aromatic side-chain atoms => 'CA' or other ring types
	'CG':  'CA',
	'CD1': 'CA',
	'CD2': 'CA',
	'CE1': 'CA',
	'CE2': 'CA',
	'CZ':  'CA',
	'CE':  'CA',
	'CH2': 'CA',
	'CZ2': 'CA',
	'CZ3': 'CA',
	'CE3': 'CA',

	# Histidine ring N
	'ND1': 'NA',
	'NE2': 'NA',

	# Tryptophan ring N
	'NE1': 'NA',

	# Some hydrogen naming
	'H2':  'H1',
	'H3':  'H1',
	'HA':  'H1',
	'HA2': 'H1',
	'HA3': 'H1',
	'HB':  'HC',
	'HB1': 'HC',
	'HB2': 'HC',
	'HB3': 'HC',
	'HD1': 'HC',
	'HD2': 'HC',
	'HD3': 'HC',
	'HE':  'HC',
	'HE1': 'HC',
	'HE2': 'HC',
	'HE3': 'HC',
	'HG':  'HC',
	'HG1': 'HC',
	'HG2': 'HC',
	'HG3': 'HC',
	'HZ':  'HC',
	'HZ1': 'HC',
	'HZ2': 'HC',
	'HZ3': 'HC',
	'HH':  'H',
	'HH11': 'H',
	'HH12': 'H',
	'HH21': 'H',
	'HH22': 'H',
}

#named tuples for atoms, residues, and interactions
from collections import namedtuple
Atom = namedtuple('Atom', [
		'id',
		'oid',
		'resi',
		'resn',
		'chain',
		'name',       # e.g. 'CZ', 'NH1', ...
		'atom_type',  # e.g. 'CT', 'N3'
		'x', 'y', 'z',
		'charge',
		'radius'
	])

	# We store a list of Atom objects plus SASA data per residue
ResidueRecord = namedtuple('ResidueRecord', [
		'chain', 'resi', 'resn',
		'atoms',                 # list of Atom
		'sasa_total',
		'sasa_polar',
		'sasa_apolar',
		'sasa_main_chain',
		'sasa_side_chain',
		'phi', 'psi', 'omega'
	])

def parse_pqr_atoms_list(atoms):
	atom_id = 0
	atoms_dict = {}
	residues = {}

	for a in atoms:
		atom_id+=1
		pqr_name = a.name.strip()
		amber_type = pqr_to_amber.get(pqr_name)
		if amber_type is None:
			# fallback guess
			if pqr_name.startswith('H'):
				amber_type = 'HC'
			elif pqr_name.startswith('O'):
				amber_type = 'O'
			elif pqr_name.startswith('N'):
				amber_type = 'N'
			elif pqr_name.startswith('C'):
				amber_type = 'CT'
			else:
				amber_type = 'CT'

		atom = Atom(
			id=atom_id,
			oid=int(a.serial),
			name=a.name,
			resn=a.res_name,
			chain=a.chain_id,
			resi=int(a.res_seq),
			atom_type=amber_type,
			x=float(a.x),
			y=float(a.y),
			z=float(a.z),
			charge=float(a.ffcharge),
			radius=float(a.radius),
		)
		atoms_dict[atom_id]=atom
		rkey = (atom.chain, atom.resi)
		if rkey not in residues:
			# create a new ResidueRecord with zero SASA, empty atom list
			new_record = ResidueRecord(
				chain           = atom.chain,
				resi            = atom.resi,
				resn            = atom.resn,
				atoms           = [atom],
				sasa_total      = 0.0,
				sasa_polar      = 0.0,
				sasa_apolar     = 0.0,
				sasa_main_chain = 0.0,
				sasa_side_chain = 0.0,
				phi             = None,
				psi             = None,
				omega           = None
			)
			residues[rkey] = new_record
		else:
			# update existing record => append this atom
			old_rec = residues[rkey]
			updated_atom_list = old_rec.atoms + [atom]
			new_rec = old_rec._replace(atoms=updated_atom_list)
			residues[rkey] = new_rec

	# Initialize bond-lists
	bonds = {}
	for a_id in atoms_dict:
		bonds[a_id] = []
	return atoms_dict, residues, bonds

		

def read_pdb_with_hydrogen(pdb_with_hydrogen):
	"""
	For each ATOM line in self.pdb_with_hydrogen, create an Atom namedtuple
	and store in self.atoms. Also update self.residues (chain, resi) -> ResidueRecord
	by appending the new atom to the 'atoms' list. By default, SASA fields = 0.
	"""
	atom_id = 0
	atoms = {}
	residues = {}
	for line in pdb_with_hydrogen:
		if line.startswith('ATOM'):
			atom_id += 1
			parts = line.split()
			if len(parts)<11:
				print(line)
			# PQR columns:
			# 0    1     2    3   4   5    6       7       8      9       10
			# ATOM oid   an  resn ch resi   x       y       z    charge   radius
			pqr_name = parts[2].strip()  # e.g. NE, CZ, ...
			# map to AMBER type
			amber_type = pqr_to_amber.get(pqr_name)
			if amber_type is None:
				# fallback guess
				if pqr_name.startswith('H'):
					amber_type = 'HC'
				elif pqr_name.startswith('O'):
					amber_type = 'O'
				elif pqr_name.startswith('N'):
					amber_type = 'N'
				elif pqr_name.startswith('C'):
					amber_type = 'CT'
				else:
					amber_type = 'CT'

			# parse residue index
			try:
				iresi = int(parts[5])
			except ValueError:
				iresi = parts[5]  # fallback if insertion code

			atom = Atom(
				id=atom_id,
				oid=int(parts[1]),
				name=pqr_name,
				resn=parts[3],
				chain=parts[4],
				resi=iresi,
				atom_type=amber_type,
				x=float(parts[6]),
				y=float(parts[7]),
				z=float(parts[8]),
				charge=float(parts[9]),
				radius=float(parts[10])
			)
			atoms[atom_id] = atom

			rkey = (atom.chain, atom.resi)
			if rkey not in residues:
				# create a new ResidueRecord with zero SASA, empty atom list
				new_record = ResidueRecord(
					chain           = atom.chain,
					resi            = atom.resi,
					resn            = atom.resn,
					atoms           = [atom],
					sasa_total      = 0.0,
					sasa_polar      = 0.0,
					sasa_apolar     = 0.0,
					sasa_main_chain = 0.0,
					sasa_side_chain = 0.0
				)
				residues[rkey] = new_record
			else:
				# update existing record => append this atom
				old_rec = residues[rkey]
				updated_atom_list = old_rec.atoms + [atom]
				new_rec = old_rec._replace(atoms=updated_atom_list)
				residues[rkey] = new_rec

	# Initialize bond-lists
	bonds = {}
	for a_id in atoms:
		bonds[a_id] = []
	return atoms, residues, bonds

def atom_dict_to_pdb(atom_dict):
	"""
	converts atoms named tuples dicts to pdb text 
	"""
	lines = []
	for atom in atom_dict.values():
		# Format: ATOM  serial  name resn chain resi   x       y       z
		line = (
			"ATOM  "
			f"{atom.oid:5d} "
			f"{atom.name:>4s} "
			f"{atom.resn:>3s} "
			f"{atom.chain:1s}"
			f"{atom.resi:4d}    "
			f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
		)
		lines.append(line)
	return "\n".join(lines)


def prepare_pdb_for_pseudo_surface(hori_atoms_dict):
	"""
	Filters the main atom dictionary to include only backbone (N, CA, C, O)
	and CB atoms (CA for GLY).
	Returns PDB-formatted lines AND an ordered list of (x,y,z) coordinates for these atoms.

	Args:
		hori_atoms_dict (dict): The Hori instance's main atoms dictionary
								(id -> Atom namedtuple).

	Returns:
		tuple: (pdb_lines_list, ordered_coords_list)
			   pdb_lines_list (list): PDB-formatted strings for backbone/CB only.
			   ordered_coords_list (list): List of (x,y,z) tuples for atoms in pdb_lines_list,
										  matching the order of atoms in the generated PDB string.
	"""
	bb_cb_atoms_list_ordered = [] # Store atom objects in a defined order

	# Ensure a consistent order of atoms, e.g., by sorting by original atom ID (oid)
	# This helps if atom_dict_to_pdb relies on the input dict's iteration order (though ideally it shouldn't)
	# or if FreeSASA processes atoms in the order they appear in the PDB file.
	sorted_atom_ids = sorted(hori_atoms_dict.keys(), key=lambda aid: hori_atoms_dict[aid].oid)

	for atom_id in sorted_atom_ids:
		atom = hori_atoms_dict[atom_id]
		atom_name = atom.name.strip().upper()
		res_name = atom.resn.strip().upper()

		is_backbone = atom_name in ["N", "CA", "C", "O"]
		is_cb = (atom_name == "CB")
		is_gly_ca_for_cb = (res_name == "GLY" and atom_name == "CA") # Treat GLY's CA as its CB

		if is_backbone or is_cb or is_gly_ca_for_cb:
			bb_cb_atoms_list_ordered.append(atom)

	if not bb_cb_atoms_list_ordered:
		print("Warning: No backbone or CB atoms found for pseudo-surface PDB preparation.")
		return [], []

	# Create a temporary dictionary for atom_dict_to_pdb if it strictly requires a dict format.
	# The keys for this temp_dict don't have to be original hori_atom_ids,
	# as long as atom_dict_to_pdb can iterate it.
	# Usingenumerate ensures a consistent order if atom_dict_to_pdb iterates based on insertion order (Python 3.7+).
	temp_dict_for_pdb_gen = {i + 1: atom for i, atom in enumerate(bb_cb_atoms_list_ordered)}

	pdb_lines_str = atom_dict_to_pdb(temp_dict_for_pdb_gen) 
	pdb_lines_list = [line + "\n" for line in pdb_lines_str.split('\n') if line.strip()]

	ordered_coords = [(atom.x, atom.y, atom.z) for atom in bb_cb_atoms_list_ordered]

	return pdb_lines_list, ordered_coords