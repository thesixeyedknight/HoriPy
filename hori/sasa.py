# hori/sasa.py
import tempfile
from .parsing import make_pdb_from_cif
from .wisz_parameters import REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED
# It's good practice to import _get_ionizable_group_type_for_sasa if it's used here,
# or ensure this logic is handled by the caller (Hori class)
# For now, to keep this function more focused on FreeSASA output,
# it will return ALL atom SASAs, and the Hori class will filter for ionizable ones.
import freesasa as fs

# fs.setVerbosity(1) # Or fs.setVerbosity(fs.nowarnings)
# It's often better to set verbosity once at the application entry point or not at all for a library.
# If you want to suppress FreeSASA warnings during normal operation:
fs.setVerbosity(1)


def compute_sasa_data(pdb_input_lines, hori_instance_residues, filetype):
	"""
	Writes PDB/CIF content to a temp file, runs freesasa,
	updates hori_instance_residues with residue-level SASA,
	and returns a dictionary of all atom-specific SASA values.

	Args:
		pdb_input_lines (list): List of strings, PDB or mmCIF atom records.
								These should represent the final structure with hydrogens.
		hori_instance_residues (dict): The Hori instance's residues dictionary.
									   This dictionary WILL BE MODIFIED IN-PLACE with residue SASA.
		filetype (str): 'pdb' or 'cif' (determines if make_pdb_from_cif is called).

	Returns:
		tuple: (all_atom_sasa_map, error_occurred_flag)
			   all_atom_sasa_map (dict): {(chain, resi_int, atom_name_str): sasa_value} for ALL atoms.
			   error_occurred_flag (bool): True if FreeSASA failed, False otherwise.
	"""
	pqr_text = "" # Renaming to reflect it's PDB-like format for FreeSASA
	if filetype == 'cif':
		# Assuming pdb_input_lines are raw CIF _atom_site lines
		pdb_formatted_lines = make_pdb_from_cif(pdb_input_lines)
		pqr_text = "".join(pdb_formatted_lines) # No need for extra newline if join handles it
	else:
		# pdb_input_lines are assumed to be PDB formatted atom lines
		pqr_text = "".join(pdb_input_lines)

	all_atom_sasa_map = {} # To store SASA for every atom parsed by FreeSASA

	if not pqr_text.strip():
		print("Error: Input to FreeSASA is empty.")
		return True # True indicates an error occurred

	with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb') as tmp:
		tmp.write(pqr_text)
		tmp.flush() # Ensure data is written before FreeSASA reads it

		try:
			structure = fs.Structure(tmp.name)
			result = fs.calc(structure)
		except RuntimeError as e: # FreeSASA can raise RuntimeError for parsing issues
			print(f"Error running FreeSASA on temporary file {tmp.name}: {e}")
			print("Temporary file content (first 1000 chars):")
			print(pqr_text[:1000] + "..." if len(pqr_text) > 1000 else pqr_text)
			return True # True indicates an error occurred
		except Exception as e: # Catch any other unexpected errors
			print(f"Unexpected error during FreeSASA calculation for {tmp.name}: {e}")
			return True


		# Process and update residue-level SASA (modifies hori_instance_residues in-place)
		res_areas = result.residueAreas()
		for chain_id_fs, chain_dict_fs in res_areas.items():
			for residue_key_str_fs, area_fs in chain_dict_fs.items():
				try:
					# freesasa residueNumber() returns string, ensure consistency
					iresi_fs = int(residue_key_str_fs)
				except ValueError:
					iresi_fs = residue_key_str_fs

				rkey_hori = (chain_id_fs.strip(), iresi_fs) # Ensure chain_id is stripped
				if rkey_hori in hori_instance_residues:
					old_rec_hori = hori_instance_residues[rkey_hori]
					updated_rec_hori = old_rec_hori._replace(
						sasa_total=area_fs.total,
						sasa_polar=area_fs.polar,
						sasa_apolar=area_fs.apolar,
						sasa_main_chain=area_fs.mainChain,
						sasa_side_chain=area_fs.sideChain
					)
					hori_instance_residues[rkey_hori] = updated_rec_hori
				# else:
				#     print(f"Warning: Residue key {rkey_hori} from FreeSASA not in Hori residues.")


		# Populate all_atom_sasa_map with SASA for every atom from the FreeSASA run
		for i in range(structure.nAtoms()):
			sasa_val = result.atomArea(i) # Corrected method call
			
			fs_atom_name = structure.atomName(i).strip()
			# fs_res_name = structure.residueName(i).strip() # For debugging
			fs_res_seq_str = structure.residueNumber(i).strip() # residueNumber() returns a string
			fs_chain_label = structure.chainLabel(i).strip()

			try:
				fs_res_seq_int = int(fs_res_seq_str)
			except ValueError:
				print(f"Warning: FreeSASA residue number '{fs_res_seq_str}' for atom {fs_atom_name} "
					  f"in chain '{fs_chain_label}' is not an integer. Skipping SASA for this atom.")
				continue
			
			atom_key = (fs_chain_label, fs_res_seq_int, fs_atom_name)
			all_atom_sasa_map[atom_key] = sasa_val
			
	return False, all_atom_sasa_map # False for no error, and the map