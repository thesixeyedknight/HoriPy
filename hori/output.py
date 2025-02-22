# output.py
import json
from datetime import datetime

def recursive_to_dict(item):
	"""
	Recursively convert namedtuples, dictionaries, lists, tuples, and sets to dicts/lists.
	"""
	if isinstance(item, tuple) and hasattr(item, '_asdict'):
		item = item._asdict()
	if isinstance(item, dict):
		return {key: recursive_to_dict(value) for key, value in item.items()}
	if isinstance(item, (list, tuple, set)):
		return [recursive_to_dict(elem) for elem in item]
	return item

def build_output_data(hori_instance):
	"""
	Build a JSON-serializable dictionary for the Hori results.
	This version omits the atomic_interactions list and stores only residue_interactions.
	"""
	metadata = {
		"fn": hori_instance.filename,       # shortened key
		"ft": hori_instance.file_type,
		"ho": hori_instance.highest_order,
		"pH": hori_instance.ph,
		"co": hori_instance.cutoff,
		"da": hori_instance.d_a_dist,
		"dha": hori_instance.dha_angle,
		"sb": hori_instance.salt_bridge_dist,
		"ppd": hori_instance.pi_pi_dist,
		"ppa": hori_instance.pi_pi_angle,
		"cpd": hori_instance.cation_pi_dist,
		"ch": list(hori_instance.chains),
		"ts": datetime.now().isoformat()
	}

	# Global atoms dictionary: id -> atom details
	atoms = {}
	for atom_id, atom in hori_instance.atoms.items():
		atoms[atom_id] = recursive_to_dict(atom)

	# Residues: store residue details and list of atom IDs
	residues = {}
	for key, residue in hori_instance.residues.items():
		key_str = f"{key[0]}_{key[1]}"
		residue_dict = recursive_to_dict(residue)
		residue_dict["aids"] = [atom["id"] for atom in residue_dict.pop("atoms", [])]
		residues[key_str] = residue_dict

	# Build residue interactions:
	# For each residue pair, process its atomic interactions.
	residue_interactions = {}
	for key, inter in hori_instance.residue_interactions.items():
		key_str = f"{key[0][0]}_{key[0][1]}__{key[1][0]}_{key[1][1]}"
		inter_dict = recursive_to_dict(inter)
		atomic_list = inter_dict.get("atomic_interactions", [])
		
		# Build a reduced list of atomic interactions, keeping only essential keys.
		reduced_atomic_list = []
		total_energy = 0.0
		for a_int in atomic_list:
			# Depending on how the atomic interaction was stored,
			# it might already be reduced (with keys "atom1_id" etc.) or still have full dicts.
			# We check and extract the id accordingly.
			if isinstance(a_int.get("atom1"), dict):
				a1_id = a_int["atom1"].get("id")
			else:
				a1_id = a_int.get("atom1_id")
			if isinstance(a_int.get("atom2"), dict):
				a2_id = a_int["atom2"].get("id")
			else:
				a2_id = a_int.get("atom2_id")
			# Create a reduced dictionary.
			reduced = {
				"a1": a1_id,
				"a2": a2_id,
				"energy": a_int.get("energy"),
				"type": a_int.get("int_type"),
				"dist": a_int.get("distance")
			}
			# Sum up the energies (if present)
			if reduced["energy"] is not None:
				total_energy += reduced["energy"]
			reduced_atomic_list.append(reduced)
		
		# Store the reduced list along with an aggregate energy value.
		inter_dict["atomic_interactions"] = reduced_atomic_list
		inter_dict["total_energy"] = total_energy
		residue_interactions[key_str] = inter_dict

	# Higher-order interactions: convert cliques (tuples) to lists.
	higher_order_interactions = {}
	for order, groups in hori_instance.higher_order_interactions.items():
		higher_order_interactions[str(order)] = [list(group) for group in groups]

	output_data = {
		"meta": metadata,
		"atoms": atoms,
		"residues": residues,
		"residue_interactions": residue_interactions,
		"higher_order_interactions": higher_order_interactions
	}
	return output_data

def save_output(hori_instance, output_filename):
	data = build_output_data(hori_instance)
	with open(output_filename, 'w') as f:
		json.dump(data, f, indent=2)
	print(f"Output successfully saved to {output_filename}")

# if __name__ == "__main__":
# 	from hori.analysis import Hori
# 	h = Hori(filename='examples/1a2p.cif', pH=13.0)
# 	save_output(h, 'analysis_output.json')
