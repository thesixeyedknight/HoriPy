"""
This file stores parameters for the Pike & Nanda (2015) empirical model
for estimating local dielectric constants.

Reference:
Pike, D. H., & Nanda, V. (2015). Empirical estimation of local dielectric
constants: Toward atomistic design of collagen mimetic peptides.
Peptide Science, 104(4), 360-370. DOI: 10.1002/bip.22644
"""
import math

# I. Atomic Polarizabilities
# --------------------------
# Extract atomic polarizability values from the Pike & Nanda paper (Table II).
# These are for AMBER united atom types (heavy atoms + bound protons).
# Units are Angstrom^3 (Å³).
# Your Hori tool's `Atom` namedtuple might use different atom type names or
# include explicit hydrogens, which are not listed here. You will need to map
# these values accordingly or decide on a strategy for explicit hydrogens.

PIKE_ATOMIC_POLARIZABILITIES = {
    # Atom Type from Paper (AMBER united) : Polarizability (Å³)
    'C':    1.382,  # General Carbon (e.g., C in GLY CA, PRO CD) - Note: Paper has C and C* with same value
    'C*':   1.382,  # Another Carbon type, same value as C in paper's table
    'C2':   1.836,  # sp2 Carbon (e.g., in aromatic rings like PHE/TYR/TRP CG,CD,CE,CZ)
    'C3':   2.222,  # sp3 Carbon (e.g., aliphatic C like ALA CB, LEU CG,CD1,CD2)
    'CA':   1.382,  # Carbonyl carbon (C in peptide bond)
    'CB':   1.529,  # Aromatic C-H C (e.g. TRP CH2, TYR CZ) - check paper context for precise mapping
    'CD':   1.768,  # Guanidinium C in ARG CZ
    'CN':   1.529,  # C in CONH2 (Asn, Gln)
    'CG':   1.768,  # C in -CH= (e.g. Pro CG) - check paper context for precise mapping
    'CH':   1.450,  # Aliphatic -CH- group C (united atom)
    'N':    1.480,  # Peptide amide Nitrogen (N)
    'N2':   1.866,  # Nitrogen in Lys NZ, Arg NE/NH1/NH2, Trp NE1
    'N3':   2.252,  # Charged Nitrogen (e.g. PRO N if charged, or specific His N) - check paper context
    'NA':   1.480,  # Nitrogen in Asn ND2, Gln NE2
    'O':    0.460,  # Carbonyl Oxygen (O in peptide bond)
    'OH':   1.050,  # Hydroxyl Oxygen (e.g., Ser OG, Thr OG1, Tyr OH)
    'O2':   0.664,  # Carboxylate Oxygen (e.g., Asp OD1/OD2, Glu OE1/OE2)
    'S':    3.2684, # Sulfur in Met SD
    'SH':   3.6643, # Sulfur in Cys SG (protonated)
    'HOH':  1.4907  # Water molecule (used for solvent contribution)
    # Note: The paper uses united atom types. If your system uses explicit hydrogens
    # (e.g., H, HC, H1, HO, HS, HN), you will need to find appropriate polarizabilities
    # from other sources or adjust the model, as this paper does not list them.
    # The examples in the original TODO (CT, N3 for Lys NZ, H, H1) need careful mapping.
    # For example, 'N3' in your example might correspond to 'N2' or 'N3' from the paper
    # depending on the specific residue and charge state context.
    # 'CT' (aliphatic C) would likely map to 'C3' or 'CH'.
}

# II. Polarizability Influence Radius
# -----------------------------------
# This is the cutoff distance (in Angstroms) beyond which the polarizability
# of a neighboring atom no longer contributes to the sum of local polarizabilities
# for the target atom.
# Found in the "Methods" section: "for all atoms within 9 Å" and
# "additive within 9 Å sphere around an atom".
POLARIZABILITY_INFLUENCE_RADIUS = 9.0  # Angstroms

# III. Parameters for Epsilon Function
# ------------------------------------
# This section is for parameters needed to convert the sum of local atomic
# polarizabilities (Σ α_local) into the local dielectric constant (ε_local).
# The paper uses Equation (2): ε = 1 + 4π N α
# Where:
#   - ε is the local dielectric constant.
#   - N is the local number density (protein atoms + water molecules) in atoms/Å³
#     within the 9 Å sphere.
#   - α (alpha_eff in our interpretation) is the sum of individual polarizabilities
#     for protein atoms AND waters in the 9 Å radius sphere.
#
# The script's structure suggests ε_local = f(Σ α_local_protein_only) or similar.
# However, the paper's 'α' in Eq. (2) already includes water contributions and
# 'N' is a local density.
#
# If we interpret the script's request as parameters for:
# ε_local = INTERCEPT + SCALE_FACTOR * (EFFECTIVE_SUMMED_POLARIZABILITY_TERM)
# where EFFECTIVE_SUMMED_POLARIZABILITY_TERM = N * α_eff
# (and α_eff = Σα_protein_atoms_in_sphere + Σα_water_molecules_in_sphere)
#
# Then the parameters are:
EPSILON_FUNCTION_PARAMS = {
    'intercept': 1.0,
    'scale_factor_for_N_alpha_product': 4.0 * math.pi,
    # Your code will need to calculate N (local density) and α_eff (total effective
    # polarizability of protein atoms + water in the 9Å sphere) for each atom.
    # The 'EFFECTIVE_SUMMED_POLARIZABILITY_TERM' would be the product (N * α_eff).
    #
    # The procedure from the paper to get N and α_eff for each target atom i:
    # 1. Σα_protein = sum of PIKE_ATOMIC_POLARIZABILITIES for protein atoms j within 9Å of atom i.
    # 2. V_sphere = (4/3)π * (POLARIZABILITY_INFLUENCE_RADIUS^3).
    # 3. V_protein_atoms_in_sphere = sum of volumes of protein atoms j (from Table II, not included here).
    # 4. V_water_space = V_sphere - V_protein_atoms_in_sphere.
    # 5. n_H2O = V_water_space / V_H2O (V_H2O from Table II, e.g., 19.1 Å³).
    # 6. α_water_total = n_H2O * PIKE_ATOMIC_POLARIZABILITIES['HOH'].
    # 7. α_eff = Σα_protein + α_water_total.
    # 8. N_atoms_in_sphere = (count of protein atoms j) + n_H2O.
    # 9. N_density = N_atoms_in_sphere / V_sphere.
    # 10. The term to use with 'scale_factor_for_N_alpha_product' is (N_density * α_eff).
}

# IV. Coulomb Constant
# --------------------
# This constant is used in Coulomb's law: E = (COULOMB_CONSTANT * q1 * q2) / (epsilon * r)
# Ensure units are consistent: if energy is kJ/mol, charges in elementary charge units (e),
# distance in Angstroms (Å), and epsilon is dimensionless.
# Your existing code uses 332.063711 * 4.184 for kJ/mol. This is standard.
# The paper uses 332 for kcal/mol.
COULOMB_CONSTANT_KJ_MOL_A_E2 = 332.063711 * 4.184 # (kcal/mol*Å/e^2) * (kJ/kcal)

# V. Other constants (if any)
# ---------------------------
# If the Pike & Nanda model uses any other specific constants for its
# dielectric calculation or application, define them here.

# The paper also uses atomic volumes for calculating N (density) and water contribution.
# These are from Table II of the paper. You'll need these if you implement the full
# N * α_eff calculation as described in EPSILON_FUNCTION_PARAMS comments.
# Example: Volume for HOH from Table II is 19.1 Å³.
# You would need the full list of atomic volumes corresponding to PIKE_ATOMIC_POLARIZABILITIES.
# For example:
PIKE_ATOMIC_VOLUMES = {
    # Atom Type from Paper : Volume (Å³) (from Table II)
    'C':    8.7,
    'C*':   8.7,
    'C2':   23.2,
    'C3':   36.7,
    'CA':   8.7,
    'CB':   8.7, # Note: Paper lists CB as 8.7, but CG/CD are larger. Check context if this is a typo or specific to united atom.
    'CD':   21.3, # e.g. ARG CZ C
    'CN':   8.7,  # e.g. ASN/GLN C in CONH2
    'CG':   21.3, # e.g. PRO CG
    'CH':   20.4, # e.g. Aliphatic CH
    'N':    13.6,
    'N2':   21.4,
    'N3':   22.7,
    'NA':   13.6,
    'O':    15.9,
    'OH':   18.0,
    'O2':   18.0,
    'S':    29.2,
    'SH':   36.7,
    'HOH':  19.1
}

# Dielectric constant of bulk solvent (epsilon_s), used in their modified GB-like term (Eq. 4).
# The paper mentions "calculated dielectric constant is 76 in a fully solvated sphere,
# consistent with the experimentally determined dielectric of bulk water." (Page 363, right col)
# Often a value like 78.5 or 80.0 is used for water at room temperature.
# The paper uses epsilon_s in Eq. 3 (referring to standard GB) and implies its use in Eq. 4 context.
BULK_SOLVENT_DIELECTRIC_EPSILON_S = 80.0 # Or a more standard 78.5 / 80.0 if preferred for general GB.

# DEFAULT_MIN_EPSILON = 1.0 # A sensible minimum for any calculated dielectric,
                            # though not explicitly a parameter from the paper's core model equations.
                            # The paper's calculated epsilons are generally >> 1.

AMBER_TO_PIKE_ATOM_TYPE_MAP = {
    # Hydrogens (P&N is united-atom, so explicit H usually not mapped for polarizability)
    'H0': None,  # H aliph. bond. to C with 1 electrwd. group (03GLY) - Implicit in P&N heavy atom
    'H': None,   # H bonded to nitrogen atoms - Implicit in P&N N types
    'HC': None,  # H aliph. bond. to C without electrwd.group - Implicit in P&N C3 or CH
    'H1': None,  # H aliph. bond. to C with 1 electrwd. group - Implicit in P&N C3 or CH
    'H2': None,  # H aliph. bond. to C with 2 electrwd.groups - Implicit in P&N C3 or CH
    'H3': None,  # H aliph. bond. to C with 3 eletrwd.groups - Implicit in P&N C3 (methyl H)
    'HA': None,  # H arom. bond. to C without elctrwd. groups - Implicit in P&N C2
    'H4': None,  # H arom. bond. to C with 1 electrwd. group - Implicit in P&N C2
    'H5': None,  # H arom. bond. to C with 2 electrwd. groups - Implicit in P&N C2
    'HO': None,  # H in hydroxyl group - Implicit in P&N OH
    'HS': None,  # H in thiol group (hydrogen bonded to sulphur) - Implicit in P&N SH
    'HW': None,  # H in TIP3P water - P&N HOH is for the whole molecule
    'HP': None,  # H bonded to C next to positively charged gr - Implicit in P&N C3 or CH

    # Carbons
    'C': 'CA',   # itp: sp2 C carbonyl group -> P&N 'CA' (Carbonyl Carbon)
    'CA': 'C2',  # itp: sp2 C pure aromatic (benzene) -> P&N 'C2' (sp2 aromatic C)
    'CB': 'C2',  # itp: sp2 aromatic C, 5&6 membered ring junction -> P&N 'C2'
    'CC': 'C2',  # itp: sp2 aromatic C, 5 memb. ring HIS -> P&N 'C2'
    'CK': 'C2',  # itp: sp2 C 5 memb.ring in purines -> P&N 'C2' (heteroaromatic)
    'CM': 'C2',  # itp: sp2 C pyrimidines in pos. 5 & 6 -> P&N 'C2' (heteroaromatic)
    'CN': 'C2',  # itp: sp2 C aromatic 5&6 memb.ring junct.(TRP) -> P&N 'C2'
                 # Note: This 'CN' is different from P&N 'CN' (C in Asn/Gln amide)
    'CQ': 'C2',  # itp: sp2 C in 5 mem.ring of purines between 2 N -> P&N 'C2' (heteroaromatic)
    'CR': 'C2',  # itp: sp2 arom as CQ but in HIS -> P&N 'C2' (heteroaromatic)
    'CT': 'C3',  # itp: sp3 aliphatic C -> P&N 'C3' (sp3 Aliphatic C). Consider P&N 'CH' if specifically methine.
    'CV': 'C2',  # itp: sp2 arom. 5 memb.ring w/1 N and 1 H (HIS) -> P&N 'C2'
    'CW': 'C2',  # itp: sp2 arom. 5 memb.ring w/1 N-H and 1 H (HIS) -> P&N 'C2'
    'C*': 'C2',  # itp: sp2 arom. 5 memb.ring w/1 subst. (TRP) -> P&N 'C2'
                 # Note: P&N also has 'C*' type, but it's "General Carbon", C2 is better for aromatic.

    # Nitrogens
    'N':  'N',   # itp: sp2 nitrogen in amide groups (peptide backbone) -> P&N 'N' (Peptide Amide N)
    'NA': 'N2',  # itp: sp2 N in 5 memb.ring w/H atom (HIS) -> P&N 'N2' (covers Trp NE1, heteroaromatic N)
    'NB': 'N2',  # itp: sp2 N in 5 memb.ring w/LP (HIS,ADE,GUA) -> P&N 'N2'
    'NC': 'N2',  # itp: sp2 N in 6 memb.ring w/LP (ADE,GUA) -> P&N 'N2'
    'N2': 'N2',  # itp: sp2 N in amino groups -> P&N 'N2' (covers Lys NZ, Arg N)
    'N3': 'N3',  # itp: sp3 N for charged amino groups (Lys, etc) -> P&N 'N3' (Charged sp3 N)
    'N*': 'N2',  # itp: sp2 N (generic) -> P&N 'N2'

    # Oxygens
    'O':  'O',   # itp: carbonyl group oxygen (peptide backbone) -> P&N 'O' (Carbonyl O)
    'OW': 'HOH', # itp: oxygen in TIP3P water -> P&N 'HOH' (Water molecule)
    'OH': 'OH',  # itp: oxygen in hydroxyl group -> P&N 'OH' (Hydroxyl O)
    'OS': 'OH',  # itp: ether and ester oxygen -> P&N 'OH' (approximate, P&N lacks specific ether/ester O)
                 # Use with caution, or map to None if a precise match is needed.
    'O2': 'O2',  # itp: carboxyl and phosphate group oxygen -> P&N 'O2' (Carboxylate O)

    # Sulfurs
    'S':  None,  # itp: S in disulfide linkage -> No direct P&N equivalent for disulfide.
                 # P&N 'S' is Met thioether, 'SH' is Cys thiol.
    'SH': 'SH',  # itp: S in cystine (assuming Cysteine thiol) -> P&N 'SH' (Thiol S)

    # Ions & Metals (No P&N parameters available)
    'Br': None,
    'C0': None,  # calcium
    'F':  None,
    'I':  None,
    'Cl': None,
    'Na': None,
    'IB': None,  # 'big ion w/ waters'
    'MG': None,
    'CU': None,
    'FE': None,
    'K':  None,
    'Rb': None,
    'Cs': None,
    'Li': None,
    'Zn': None,

    # Phosphorus (No P&N parameters available)
    'P':  None,

    # Other Water Model Atoms
    'OW_spc':     'HOH', # SPC Water OW -> P&N 'HOH'
    'HW_spc':     None,  # SPC Water HW - Implicit in P&N HOH
    'OW_tip4pew': 'HOH', # tip4pEW OW -> P&N 'HOH'
    'HW_tip4pew': None,  # tip4pEW HW - Implicit in P&N HOH
    'OW_tip4p':   'HOH', # tip4p OW -> P&N 'HOH'
    'HW_tip4p':   None,  # tip4p HW - Implicit in P&N HOH
    'OW_tip5p':   'HOH', # tip5p OW -> P&N 'HOH'
    'HW_tip5p':   None,  # tip5p HW - Implicit in P&N HOH

    # Virtual/Dummy Sites (No polarizability)
    'MW':   None,
    'MCH3': None,
    'MNH3': None,
}

# Bundle them for easier import/use
PIKE_NANDA_ALL_PARAMS = {
    "PIKE_ATOMIC_POLARIZABILITIES": PIKE_ATOMIC_POLARIZABILITIES,
    "PIKE_ATOMIC_VOLUMES": PIKE_ATOMIC_VOLUMES,
    "POLARIZABILITY_INFLUENCE_RADIUS": POLARIZABILITY_INFLUENCE_RADIUS,
    "EPSILON_FUNCTION_PARAMS": EPSILON_FUNCTION_PARAMS,
    "COULOMB_CONSTANT_KJ_MOL_A_E2": COULOMB_CONSTANT_KJ_MOL_A_E2,
    "AMBER_TO_PIKE_ATOM_TYPE_MAP": AMBER_TO_PIKE_ATOM_TYPE_MAP,
    "BULK_SOLVENT_DIELECTRIC_EPSILON_S": BULK_SOLVENT_DIELECTRIC_EPSILON_S
}