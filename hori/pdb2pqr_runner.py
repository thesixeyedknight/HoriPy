import tempfile
import subprocess as sp
import argparse
from mod_pdb2pqr.pdb import ATOM
from mod_pdb2pqr.main import main_driver

def run_pdb2pqr(atom_lines, file_type, ph=7.0):
	"""
	Runs pdb2pqr on a list of atom lines and returns the PQR output.

	Parameters:
	- atom_lines (list of str): Lines from a PDB or mmCIF file starting with ATOM/HETATM.
	- file_type (str): Type of the file ('pdb' or 'cif').

	Returns:
	- pqr_content (str): The generated PQR content.
	- missed_residues (list): List of residues that were missed or had issues.
	- pka_df (pandas.DataFrame or None): DataFrame containing pKa information, if any.
	- biomolecule (Biomolecule): The biomolecule object processed by pdb2pqr.
	"""
	# Step 1: Convert atom lines to pdblist using UnifiedATOM
	pdblist_for_main = []
	for line in atom_lines:
		pdblist_for_main.append(ATOM(line=line, filetype=file_type))
	print(f'Passing {len(pdblist_for_main)} atoms to main_driver')
	# Step 2: Create an argparse.Namespace with required arguments
	args = argparse.Namespace(
		input_path='custom_input',       # Dummy value since we're using a custom PDB list
		output_pqr='dummy_output.pqr',   # Dummy value; output won't be used if modified
		log_level='ERROR',            # Suppress logging
		ff='AMBER',                       # Default forcefield; adjust as needed
		userff=None,
		clean=False,
		debump=True,
		opt=True,
		keep_chain=True,
		assign_only=False,
		ffout=None,
		usernames=None,
		apbs_input=None,
		pdb_output=None,
		ligand=None,
		whitespace=False,
		neutraln=False,
		neutralc=False,
		drop_water=False,
		include_header=False,
		pka_method=None,
		ph=ph,
		parameters=None,                  # Default parameters; adjust if necessary
	)
	# Step 3: Call main_driver with the custom PDB list
	results, bio_mol = main_driver(args, custom_pdblist=pdblist_for_main)
	matched_atoms = results.get('matched_atoms')
	
	return matched_atoms








###########################################################

def run_pdb2pqr_old(pdb_for_haad):
	pdb_text = "".join(pdb_for_haad)
	with tempfile.NamedTemporaryFile(mode='w', suffix='.pdb') as tmp_in, \
			tempfile.NamedTemporaryFile(mode='r', suffix='.pqr') as tmp_out:
		tmp_in.write(pdb_text)
		tmp_in.flush()

		cmd = ['pdb2pqr', '--ff=AMBER', '--keep-chain', tmp_in.name, tmp_out.name]
		try:
			proc = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True)
			_, stderr = proc.communicate()
			if proc.returncode != 0:
				print("pdb2pqr error:", stderr)
				return
			tmp_out.seek(0)
			return tmp_out.read().splitlines()
		except Exception as e:
			print("Exception in pdb2pqr:", e)

