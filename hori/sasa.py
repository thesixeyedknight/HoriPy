import tempfile
from .parsing import make_pdb_from_cif
import freesasa as fs
fs.setVerbosity(1)

def compute_residue_sasa(pdb_with_hydrogen, residues, filetype):
	"""
	Writes the entire PQR to a temp file, runs freesasa, then merges
	SASA data back into self.residues. We match chain & residue number
	with FreedSASA output, ignoring leftover or new residues if any.
	"""
	if filetype == 'cif':
		lines = make_pdb_from_cif(pdb_with_hydrogen)
		pqr_text = "".join(lines) + "\n"
	else:
		pqr_text = "".join(pdb_with_hydrogen)
	with tempfile.NamedTemporaryFile(mode='w', suffix='.cif') as tmp:
		tmp.write(pqr_text)
		tmp.flush()

		structure = fs.Structure(tmp.name)
		result = fs.calc(structure)
		res_areas = result.residueAreas()
		# res_areas = { chain_id: { residue_key_str: ResidueArea, ... }, ... }

	# Merge SASA data
	for chain_id, chain_dict in res_areas.items():
		for residue_key_str, area in chain_dict.items():
			# FreedSASA residue_key is typically a string, e.g. "317"
			# Attempt to parse it
			try:
				iresi = int(residue_key_str)
			except ValueError:
				iresi = residue_key_str

			rkey = (chain_id, iresi)
			if rkey in residues:
				old_rec = residues[rkey]
				# keep same atoms & resn, just update SASA
				updated_rec = old_rec._replace(
					sasa_total      = area.total,
					sasa_polar      = area.polar,
					sasa_apolar     = area.apolar,
					sasa_main_chain = area.mainChain,
					sasa_side_chain = area.sideChain
				)
				residues[rkey] = updated_rec
