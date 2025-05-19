# hori/wisz_parameters.py

"""
This file stores the parameters for the Wisz and Hellinga (2003)
empirical model for electrostatic interactions.
"""

# Parameters from Table II of Wisz & Hellinga (2003)
# PROTEINS: Structure, Function, and Genetics 51:360-377 (2003)
WISZ_HELLINGA_PARAMS = {
	# Solvation parameters
	"solv_epsilon_C": 19.0,  # Dielectric for core desolvation
	"solv_epsilon_B": 32.0,  # Dielectric for boundary desolvation
	"R_Born": 2.0,           # Born radius (Angstroms)

	# Pairwise charge-charge (QQ) interaction parameters
	# QQ_alpha: dielectric constant at large separation
	# QQ_lambda: damping constant for distance-dependent dielectric
	"QQ_alpha_C": 110.0,
	"QQ_alpha_B": 110.0,
	"QQ_alpha_S": 150.0,
	"QQ_lambda_C": 0.2,
	"QQ_lambda_B": 0.4,
	"QQ_lambda_S": 1.4,

	# Pairwise charge-polar (QP) interaction parameters
	# QP_alpha: dielectric constant at large separation
	# QP_lambda: damping constant for distance-dependent dielectric
	"QP_alpha_C": 90.0,
	"QP_alpha_B": 110.0,
	"QP_alpha_S": 150.0,
	"QP_lambda_C": 0.2,
	"QP_lambda_B": 0.5,
	"QP_lambda_S": 0.9,

	# Charge-polar (HbQP) hydrogen bonding parameters (kcal/mol)
	"Delta_HbQP_G_C": -1.0,
	"Delta_HbQP_G_B": -0.25,
	"Delta_HbQP_G_S": 1.0,

	# Charge-charge (HbQQ) hydrogen bonding parameters (kcal/mol)
	"Delta_HbQQ_G_C": -5.25,
	"Delta_HbQQ_G_B": -2.25,
	"Delta_HbQQ_G_S": 3.75,

	# Ion pairing parameters (dielectric constants)
	"IP_epsilon_C": 16.0,
	"IP_epsilon_B": 135.0,
	"IP_epsilon_S": 165.0,
}

# Reference SASA values for ionizable groups in extended tripeptides (GXG)
# User needs to ensure these values are accurate and comprehensive.
# Atoms listed are those considered "ionizable atoms" for SASA calculation for that group.
REFERENCE_IONIZABLE_GROUP_SASA_EXTENDED = {
	"ASP_SIDECHAIN": {"atoms": ["OD1", "OD2"], "sasa": 78.21},
	"GLU_SIDECHAIN": {"atoms": ["OE1", "OE2"], "sasa": 88.7},
	"LYS_SIDECHAIN": {"atoms": ["NZ"], "sasa": 57.03},
	"ARG_SIDECHAIN": {"atoms": ["NE", "NH1", "NH2", "CZ"], "sasa": 131.66},
	"HIS_SIDECHAIN": {"atoms": ["ND1", "NE2", "CG", "CD2", "CE1"], "sasa": 124.54}, # Note: Paper might classify ND1/NE2 differently. User specified these atoms.
	"CYS_SIDECHAIN": {"atoms": ["SG"], "sasa": 61.62}, # For reduced Cys
	"TYR_SIDECHAIN": {"atoms": ["OH", "CZ"], "sasa": 51.45}, # Paper uses OH for Tyr ionization. CZ included per user's list.
	"NTERM": {"atoms": ["N"], "sasa": 50.31}, # Standard N-terminus
	"NTERM_PRO": {"atoms": ["N"], "sasa": 22.07}, # N-terminus if Proline
	"CTERM": {"atoms": ["O", "OXT"], "sasa": 43.44}, # Standard C-terminus
	# Heme propionate is mentioned in the paper but not standard.
	# If needed, it should be added here. Example:
	# "HEME_PROP_A": {"atoms": ["O1A", "O2A"], "sasa": VALUE},
	# "HEME_PROP_D": {"atoms": ["O1D", "O2D"], "sasa": VALUE},
}

# Default values for calculations if not specified elsewhere
DEFAULT_WATER_DIELECTRIC = 80.0 # Epsilon_w
DEFAULT_IONIC_STRENGTH = 0.0 # Molar, if not specified by experiment
DEFAULT_TEMPERATURE = 298.15 # Kelvin (25 C)
# Physical constants (ensure units are consistent with Wisz paper: kcal/mol, Angstroms, electron charge units)
# Conversion factor 332 used in paper implies charges are in electron units, distance in Angstroms, energy in kcal/mol
# k_e = 332.0637 (kcal/mol * A / e^2) if epsilon is unitless.
# The paper uses 332.
# R (Boltzmann constant in kcal/mol·K) = 0.0019872041
# RT at 298.15 K = 0.0019872041 * 298.15 = 0.592 kcal/mol
# 2.303 * RT = 1.364 kcal/mol at 298.15 K (used for pKa calculations)

# Born radii (R_Born) for different groups are mentioned as a single value in Table II
# but in some models they can be group-specific. The paper uses one R_Born.

# Partial charges for non-titrating groups (PARSE/CHARMM22)
# The paper mentions using PARSE for neutral amino acids and CHARMM22 for non-AA groups.
# For HoriPy, these are already handled by the force field selected for the PQR generation (AMBER in current HoriPy).
# The Wisz model itself would need these partial charges for the QP (charge-polar) term.
#HoriPy's `Atom` object (derived from PDB2PQR output) already has reliable partial charges, we can use those.

# Charges of ionizable groups in their neutral and charged states:
# The paper states: "partial charges for the neutral form of these side chains are derived by distributing a unit charge
# equally over the two carboxylate oxygens [for Asp/Glu/C-term], while leaving the group unprotonated in the calculations."
# "For basic groups... distributing the subtraction of a unit charge equally over the three hydrogens on each of these groups,
# while leaving the site fully protonated in the calculations."
# "The charges of all ionizable groups except histidine are localized to one or two atoms of the side chain,
# obtained by subtracting or adding one unit of charge from the neutral form of an acidic or basic residue"
# This implies that for the Q_i Q_j calculations involving ionizable groups, q_i and q_j are +1 or -1 (or fractional if pH-dependent).
# For Q_i q_j (charge-polar), q_j are the permanent partial charges.