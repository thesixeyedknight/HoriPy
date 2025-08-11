from .math_utils import distance
from multiprocessing import Pool
from collections import namedtuple
from .classify_interactions import classify_interaction
from .energy import compute_interaction_energy

def populate_bonds(residues, bonds):
	"""
	Build the bond map for each residue using the residue_hydrogen_mapping and
	residue_heavy_atom_bonds dictionaries for known residue types. Each hydrogen
	or heavy atom bond is added if found in the respective dictionaries.

	If a residue type is not in the dictionary, fallback logic is applied:
	- Hydrogens: Nearest heavy atom within 1.5 Å.
	- Heavy atoms: Nearest heavy atom within 2 Å for unlisted bonds.
	"""

	# Dictionary to map hydrogens to heavy atoms in the same residue
	residue_hydrogen_mapping = {
		"ALA": {"H": "N", "HA": "CA", "HB1": "CB", "HB2": "CB", "HB3": "CB"},
		"ARG": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG", "HD2": "CD", "HD3": "CD", "HE": "NE", "HH11": "NH1", "HH12": "NH1", "HH21": "NH2", "HH22": "NH2"},
		"ASN": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD21": "ND2", "HD22": "ND2"},
		"ASP": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB"},
		"CYS": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG": "SG"},
		"GLN": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG", "HE21": "NE2", "HE22": "NE2"},
		"GLU": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG"},
		"GLY": {"H": "N", "HA2": "CA", "HA3": "CA"},
		"HIS": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "ND1", "HE1": "CE1", "HD2": "CD2", "HE2": "CE2"},
		"ILE": {"H": "N", "HA": "CA", "HB": "CB", "HG12": "CG1", "HG13": "CG1", "HG2": "CG2", "HD11": "CD1", "HD12": "CD1", "HD13": "CD1"},
		"LEU": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG": "CG", "HD1": "CD1", "HD2": "CD2"},
		"LYS": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG", "HD2": "CD", "HD3": "CD", "HE2": "CE", "HE3": "CE", "HZ1": "NZ", "HZ2": "NZ", "HZ3": "NZ"},
		"MET": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG", "HE1": "CE", "HE2": "CE", "HE3": "CE"},
		"PHE": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HD2": "CD2", "HE1": "CE1", "HE2": "CE2", "HZ": "CZ"},
		"PRO": {"H2": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG2": "CG", "HG3": "CG", "HD2": "CD", "HD3": "CD"},
		"SER": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HG": "OG"},
		"THR": {"H": "N", "HA": "CA", "HB": "CB", "HG1": "OG1", "HG21": "CG2", "HG22": "CG2", "HG23": "CG2"},
		"TRP": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HE1": "NE1", "HE3": "CE3", "HZ2": "CZ2", "HZ3": "CH2", "HH2": "NH2"},
		"TYR": {"H": "N", "HA": "CA", "HB2": "CB", "HB3": "CB", "HD1": "CD1", "HD2": "CD2", "HE1": "CE1", "HE2": "CE2", "HH": "OH"},
		"VAL": {"H": "N", "HA": "CA", "HB": "CB", "HG11": "CG1", "HG12": "CG1", "HG13": "CG1", "HG21": "CG2", "HG22": "CG2", "HG23": "CG2"},
	}

	# Dictionary to define heavy atom bonds for each residue
	residue_heavy_atom_bonds = {
		"ALA": [("N", "CA"), ("CA", "CB"), ("CA", "C"), ("C", "O")],
		"ARG": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")],
		"ASN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
		"ASP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
		"CYS": [("N", "CA"), ("CA", "CB"), ("CB", "SG")],
		"GLN": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")],
		"GLU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
		"GLY": [("N", "CA"), ("CA", "C"), ("C", "O")],
		"HIS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("ND1", "CE1"), ("CE1", "NE2"), ("NE2", "CD2"), ("CD2", "CG")],
		"ILE": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1"), ("CA", "C"), ("C", "O")],
		"LEU": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CA", "C"), ("C", "O")],
		"LYS": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ"),("CA", "C"), ("C", "O")],
		"MET": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE"), ("CA", "C"), ("C", "O")],
		"PHE": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"),("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CA", "C"), ("C", "O")],
		"PRO": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "N"), ("CA", "C"), ("C", "O")],
		"SER": [("N", "CA"), ("CA", "CB"), ("CB", "OG"), ("CA", "C"), ("C", "O")],
		"THR": [("N", "CA"), ("CA", "CB"), ("CB", "OG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")],
		"TRP": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CD1", "NE1"), ("NE1", "CE2"),("CE2", "CD2"), ("CD2", "CG"), ("CE2", "CZ2"), ("CZ2", "CH2"), ("CH2", "CZ3"),("CZ3", "CE3"), ("CE3", "CD2"), ("CA", "C"), ("C", "O")],
		"TYR": [("N", "CA"), ("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"), ("CD1", "CE1"),("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH"), ("CA", "C"), ("C", "O")],
		"VAL": [("N", "CA"), ("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CA", "C"), ("C", "O")]
	}

	for (chain, resi), res_rec in residues.items():
		resn = res_rec.resn.upper()
		atoms_in_res = res_rec.atoms

		if resn in residue_hydrogen_mapping:
			# 1) Handle hydrogen bonds using the dictionary-based approach
			h_map = residue_hydrogen_mapping[resn]

			for atom in atoms_in_res:
				if atom.name in h_map:
					heavy_name = h_map[atom.name]
					# Find the heavy atom in the same residue
					heavy_atom = next((a for a in atoms_in_res if a.name == heavy_name), None)
					if heavy_atom is not None:
						# Bond them, avoiding redundancies
						if heavy_atom.id not in bonds[atom.id]:
							bonds[atom.id].append(heavy_atom.id)
						if atom.id not in bonds[heavy_atom.id]:
							bonds[heavy_atom.id].append(atom.id)

		else:
			# Fallback: Nearest heavy atom for hydrogens
			heavy_atoms = [a for a in atoms_in_res if not a.atom_type.startswith('H')]
			hydrogens = [a for a in atoms_in_res if a.atom_type.startswith('H')]

			for h in hydrogens:
				best_dist = 999.0
				best_heavy = None
				for hvy in heavy_atoms:
					dist = distance(h, hvy)
					if dist < best_dist:
						best_dist = dist
						best_heavy = hvy
				# If below threshold, bond them
				if best_heavy and best_dist < 1.5:
					if best_heavy.id not in bonds[h.id]:
						bonds[h.id].append(best_heavy.id)
					if h.id not in bonds[best_heavy.id]:
						bonds[best_heavy.id].append(h.id)

		if resn in residue_heavy_atom_bonds:
			# 2) Handle heavy atom bonds using the dictionary-based approach
			heavy_bonds = residue_heavy_atom_bonds[resn]

			for bond_pair in heavy_bonds:
				atom1 = next((a for a in atoms_in_res if a.name == bond_pair[0]), None)
				atom2 = next((a for a in atoms_in_res if a.name == bond_pair[1]), None)

				if atom1 and atom2:
					# Bond them, avoiding redundancies
					if atom2.id not in bonds[atom1.id]:
						bonds[atom1.id].append(atom2.id)
					if atom1.id not in bonds[atom2.id]:
						bonds[atom2.id].append(atom1.id)

		else:
			# Fallback: Nearest heavy atom bonds for heavy atoms not in dictionary
			heavy_atoms = [a for a in atoms_in_res if not a.atom_type.startswith('H')]

			for i, atom1 in enumerate(heavy_atoms):
				for atom2 in heavy_atoms[i + 1:]:
					dist = distance(atom1, atom2)
					if dist < 2.0:
						if atom2.id not in bonds[atom1.id]:
							bonds[atom1.id].append(atom2.id)
						if atom1.id not in bonds[atom2.id]:
							bonds[atom2.id].append(atom1.id)

def compute_distance_residue_pairs(args):
	"""
	Compute distances between all atoms in two residues and return distances <= cutoff.
	"""
	res_rec1, res_rec2, cutoff, distance_fn = args
	result = {}
	for a1 in res_rec1.atoms:
		for a2 in res_rec2.atoms:
			dist = distance_fn(a1, a2)
			if dist <= cutoff:
				key = (min(a1.id, a2.id), max(a1.id, a2.id))
				result[key] = dist
	return result

def build_distance_map_parallel(residues, distance_map, cutoff=7.0, num_cores=1):
	"""
	Precompute distances for all pairs of atoms across residues using parallel processing.
	"""
	# Group residues
	all_residues_list = list(residues.values())
	
	# Create a list of residue pairs to process
	residue_pairs = []
	for i, res_rec1 in enumerate(all_residues_list):
		for res_rec2 in all_residues_list[i+1:]:
			if res_rec1.chain == res_rec2.chain and abs(res_rec1.resi - res_rec2.resi) < 2:
				continue
			residue_pairs.append((res_rec1, res_rec2, cutoff, distance))
	
	# Use multiprocessing Pool to parallelize computation
	with Pool(processes=num_cores) as pool:
		results = pool.map(compute_distance_residue_pairs, residue_pairs)
	
	# Combine results into the distance_map
	for result in results:
		distance_map.update(result)
	print(f"Distance map size: {len(distance_map)}, cutoff={cutoff}")

Interaction = namedtuple('Interaction', [
		'atom1', 'atom2',
		'distance',
		'int_type',
		'energy',
		'kbp_energy',
		'geom_metrics'
	])

def find_atomic_interactions(distance_map, atoms, bonds, residues, atom_interactions, amber_nonbonded, user_params, hori_instance=None, kbp_manager=None):
	for (id1, id2), dist in distance_map.items():
		a1 = atoms[id1]
		a2 = atoms[id2]

		interaction_details = classify_interaction(distance_map, a1, a2, dist, atoms, bonds, residues, user_params)
		if interaction_details:
			itype = interaction_details.pop('type')
			geom_metrics = interaction_details
			en, kbp_en = compute_interaction_energy(a1, a2, dist, itype, residues, atoms, bonds, amber_nonbonded, user_params, hori_instance=hori_instance, kbp_manager=kbp_manager, geom_metrics=geom_metrics)
			inter = Interaction(a1, a2, dist, itype, en, kbp_en, geom_metrics)
			key = (min(a1.id, a2.id), max(a1.id, a2.id))
			if itype == 'salt_bridge':
				if en > -2.0: #Ensure that salt bridges are above threshold energy
					continue
			atom_interactions[key] = inter

def find_residue_interactions(atom_interactions, residue_interactions):
	"""
	From atomic interactions, find residue-residue pairwise interactions.
	Store residues and corresponding atomic interactions.
	"""
	for key, interaction in atom_interactions.items():
		res1 = (interaction.atom1.chain, interaction.atom1.resi)
		res2 = (interaction.atom2.chain, interaction.atom2.resi)
		if res1 == res2:
			continue  # Ignore intra-residue interactions
			
		# Maintain consistency in key ordering
		res_pair = tuple(sorted([res1, res2]))
		
		if res_pair not in residue_interactions:
			residue_interactions[res_pair] = {
				"atomic_interactions": [],
				"int_types": set(),
			}
		
		residue_interactions[res_pair]["atomic_interactions"].append(key)
		residue_interactions[res_pair]["int_types"].add(interaction.int_type)