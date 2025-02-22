from pathlib import Path

DATA_DIR = Path(__file__).parent.parent/"data"
ATOMTYPES_FILE = DATA_DIR / "atomtypes.atp"
FFNONBONDED_FILE = DATA_DIR / "ffnonbonded.itp"

def parse_amber_atomtypes(amber_atomtypes, atomtypes_path=ATOMTYPES_FILE):
	try:
		with open(atomtypes_path, 'r') as f:
			for line in f:
				line = line.strip()
				if not line or line.startswith(';') or line.startswith('#'):
					continue
				parts = line.split()
				atomtype = parts[0]
				mass = float(parts[1])
				comment = None
				if ';' in line:
					comment = line.split(';', 1)[1].strip()
				amber_atomtypes[atomtype] = (mass, comment)
	except IOError:
		print(f"Error reading {atomtypes_path}")

def parse_amber_ffnonbonded(amber_nonbonded, ffnonbonded_path=FFNONBONDED_FILE):
	in_atomtypes = False
	try:
		with open(ffnonbonded_path, 'r') as f:
			for line in f:
				line = line.strip()
				if not line or line.startswith(';') or line.startswith('#'):
					continue
				if line.lower().startswith('[ atomtypes ]'):
					in_atomtypes = True
					continue
				if in_atomtypes:
					if line.startswith('['):
						# end of block
						break
					parts = line.split()
					if len(parts) < 7:
						continue
					name = parts[0]
					sigma = float(parts[5])
					epsilon = float(parts[6])
					amber_nonbonded[name] = (sigma, epsilon)
	except IOError:
		print(f"Error reading {ffnonbonded_path}")