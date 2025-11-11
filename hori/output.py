# output.py
import json
from datetime import datetime
import numpy as np
import dataclasses

def recursive_to_dict(item):
	"""
	Recursively convert namedtuples, dictionaries, lists, tuples, and sets to dicts/lists. Includes numpy types.
	"""
	# Check for numpy types FIRST
	if isinstance(item, (np.floating, np.float64)):
		return float(item)
	if isinstance(item, (np.integer, np.int64)):
		return int(item)
	if isinstance(item, np.ndarray):
		return item.tolist()
	if dataclasses.is_dataclass(item):
		item = dataclasses.asdict(item)
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
	This version stores atomic interactions at the top level for data integrity
	and uses keys for referencing within residue interactions.
	"""
	metadata = {
		"version": "5.0",
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

	# build atomic interactions:
	atomic_interactions = {}
	for key, interaction in hori_instance.atom_interactions.items():
		key_str = f"{key[0]}_{key[1]}"
		atomic_interactions[key_str] = {
			"a1": interaction.atom1.id,
			"a2": interaction.atom2.id,
			"energy": interaction.energy,
			"kbp_energy": interaction.kbp_energy,
			"type": interaction.int_type,
			"dist": interaction.distance,
			"geom_metrics": recursive_to_dict(interaction.geom_metrics),
			"nis": recursive_to_dict(interaction.nis)
		}
	# Build residue interactions:
	residue_interactions = {}
	for res_pair_key, res_inter_data in hori_instance.residue_interactions.items():
		key_str = f"{res_pair_key[0][0]}_{res_pair_key[0][1]}__{res_pair_key[1][0]}_{res_pair_key[1][1]}"
		atomic_keys = res_inter_data["atomic_interactions"]
		total_energy = sum(hori_instance.atom_interactions[key].energy or 0 for key in atomic_keys)
		residue_interactions[key_str] = {
			"atomic_interaction_keys": [f"{k[0]}_{k[1]}" for k in atomic_keys],
			"int_types": list(res_inter_data["int_types"]),
			"total_energy": total_energy
		}

	# Higher-order interactions: convert cliques (tuples) to lists.
	higher_order_interactions = {}
	for order, groups in hori_instance.higher_order_interactions.items():
		higher_order_interactions[str(order)] = [list(group) for group in groups]

	output_data = {
		"meta": metadata,
		"atoms": atoms,
		"residues": residues,
		"atomic_interactions": atomic_interactions,
		"residue_interactions": residue_interactions,
		"higher_order_interactions": higher_order_interactions
	}
	return output_data

def save_output(hori_instance, output_filename):
	data = build_output_data(hori_instance)
	with open(output_filename, 'w') as f:
		json.dump(data, f, indent=2)
	print(f"Output successfully saved to {output_filename}")
