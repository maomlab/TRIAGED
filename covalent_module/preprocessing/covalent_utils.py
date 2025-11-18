import os
import random
import string
import pickle
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from pdbeccdutils.core.component import ConformerType 

WARHEAD_REACTIONS = { "nitrile": "[C:3][C:4]#[N:5]>>[C:3][13C:4]=[N:5]", 
"alkylhalide" : "[CX4;CH,CH2:2][I,Br,Cl:3]>>[13C:2]",
"vinyl-sulfone" : "[C:3]=[C:4][S:5](=O)=O>>[13C:3][C:4][S:5](=O)=O", # for CYS rxn; might be diff for HIS (Schneider, Grabowsky 2015)
"acrylamide" : "[C:2]=[C:3]-C(=O)-[N:4]>>[13C:3]-[C:2]-C(=O)-[N:4]",
"nitrile2": "[N:4]#[C:5]>>[N:4]=[13C:5]"
}
WARHEAD_REACTANTS = {name: smarts.split(">>")[0] for name, smarts in WARHEAD_REACTIONS.items()}
compiled_reactants = {name: Chem.MolFromSmarts(smarts) for name, smarts in WARHEAD_REACTANTS.items()}

# use ccd_pkl conda env
def verify_covalent(res_name):
    '''
    Checks if a residue is a covalent amino acid.
    :param res_name: Three-letter residue code (e.g., 'CYS').
    :return (bool): True if covalent-capable, False otherwise.
    '''
    cov_aa = {"CYS", "SER", "THR", "LYS", "HIS", "PRO", "TYR", "GLU", "ASP"}
    return res_name.upper() in cov_aa

def residue_cov_atom(res_name):
    '''
    Given a residue name, return the atom commonly used for covalent bonding.
    :param res_name: Three-letter residue code.
    :return: Atom name string.
    '''
    mapping = {
        "CYS": "SG",
        "SER": "OG",
        "THR": "OG1",
        "LYS": "NZ",
        "HIS": "NE2",  # could also be NE1, ND1
        "PRO": "N",
        "TYR": "OH",
        "GLU": "OE2",
        "ASP": "OD2",
    }
    try:
        return mapping[res_name]
    except KeyError:
        raise ValueError(f"[ERROR] Residue '{res_name}' not supported for covalent docking.")

def identify_warhead(smiles):
    '''
    Identifies and returns the names of covalent warheads found in a given molecule.
    :param smiles: str
        The SMILES representation of the molecule.
    :returns: list
        A list of warhead names found in the molecule.
    '''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("[ERROR] Invalid SMILES string.")
        return []
        
    found_warheads = []
    for name, pattern in compiled_reactants.items():
        if pattern and mol.HasSubstructMatch(pattern):
            found_warheads.append(name)

    if len(found_warheads) > 1:
        if VERBOSE: print("[WARNING] More than 1 warhead found. Choosing first match.")
    elif len(found_warheads) == 0:
        raise ValueError(f'[ERROR] No matching warhead was found for {smiles}')
    return found_warheads[0]

def ligand_cov_atom(no_lg_smiles):
    '''
    Finds the index of the first Carbon-13 ([13C]) atom in a SMILES string.
    :param smiles_string: str
        The SMILES string to parse.
    :returns: int
        The 0-based index of the first C13 atom found.
        Returns -1 if no C13 atom is found or the SMILES is invalid.
    '''
    mol = Chem.MolFromSmiles(no_lg_smiles)

    if not mol:
        print("[ERROR] Invalid SMILES string provided.")
        return -1

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetIsotope() == 13:
            return f'C{atom.GetIdx()}'
    return -1

def remove_leaving_group(smiles):
    '''
    Apply a covalent warhead-specific reaction to remove the leaving group
    from a ligand SMILES string.
    :param smiles: str
        Input ligand SMILES string.

    :returns: tuple
            smi_no_c13: str
                Processed SMILES with leaving group removed and isotopic labels stripped.
            lig_atom: int
                Index/identifier of the covalent attachment atom.
            wh_type: str
                Identified warhead type. 
    '''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"

    mol = Chem.MolFromSmiles(smiles)
    wh_type = identify_warhead(smiles)

    rxn = AllChem.ReactionFromSmarts(WARHEAD_REACTIONS[wh_type])

    products = rxn.RunReactants((mol,))

    for prod_tuple in products: # returns ((<rdkit.Chem.rdchem.Mol object at 0x14d267c24040>,),) 
        for prod in prod_tuple:
            smi_no_lg = Chem.MolToSmiles(prod)

    lig_atom = ligand_cov_atom(smi_no_lg)
    smi_no_c13 = smi_no_lg.replace("13", "")

    if VERBOSE: print('[SUCCESS] Leaving group removed:', smi_no_c13)

    return smi_no_c13, lig_atom, wh_type

def compute_3d(mol) -> bool:
    '''Generate 3D coordinates using EKTDG method.
    Taken from `pdbeccdutils.core.component.Component`.
    :param mol: Mol
        The RDKit molecule to process
    :return: bool
        Whether computation was successful.
    '''
    options = rdkit.Chem.AllChem.ETKDGv3()
    options.clearConfs = False
    conf_id = -1

    try:
        conf_id = rdkit.Chem.AllChem.EmbedMolecule(mol, options)
        rdkit.Chem.AllChem.UFFOptimizeMolecule(mol, confId=conf_id, maxIters=1000)

    except RuntimeError:
        pass  # Force field issue here
    except ValueError:
        pass  # sanitization issue here

    if conf_id != -1:
        conformer = mol.GetConformer(conf_id)
        conformer.SetProp("name", ConformerType.Computed.name)
        conformer.SetProp("coord_generation", f"ETKDGv3")
        
        return True
    return False

def get_conformer(mol, c_type: ConformerType):
    ''' 
    Retrieve an rdkit object for a deemed conformer.
    Taken from `pdbeccdutils.core.component.Component`.
    :param mol: Mol
        The molecule to process.
    :param c_type: ConformerType
        The conformer type to extract.

    :return: Conformer
        The desired conformer, if any.
    '''
    for c in mol.GetConformers():
        try:
            if c.GetProp("name") == c_type.name:
                return c
        except KeyError:  # noqa: PERF203
            pass

    msg = f"Conformer {c_type.name} does not exist."
    raise ValueError("[ERROR] ", msg)

def process_covalent_smiles(ccd_db, smiles, compound_id):
    '''
    Creates a covalent CCD ligand from a conanical SMILES string using Rdkit.
    :param smiles: str
        Sanitized canonical SMILES string of the covalent ligand.
    :param compound_id: str
        Unique identifier for the compound, used in naming the CCD file.
        If not given, a random CCD code will be generated.
    :param ccd_db: str
        Path to the CCD database directory where .pkl files are stored.

    :return: str
        Unique CCD code (or the vault_id) for the covalent ligand pkl file.
    '''
    VERBOSE = os.environ.get("VERBOSE", "FALSE").upper() == "TRUE"
    DEBUG = os.environ.get("DEBUG", "FALSE").upper() == "TRUE"

    pkl_file = os.path.join(ccd_db, f"{compound_id}.pkl")

    if DEBUG and os.path.exists(pkl_file): 
        if VERBOSE: print("[WARNING] {pkl_file} will be removed and will be re-written!")
        os.remove(pkl_file)
        
    if not os.path.exists(pkl_file):
        mol_sdf = Chem.MolFromSmiles(smiles) 
        success = compute_3d(mol_sdf)
        if not success:
            print(f"[ERROR] 3D coordinate generation failed for the provided SMILES: {smiles}.")
        _ = get_conformer(mol_sdf, ConformerType.Computed)
        for _, atom in enumerate(mol_sdf.GetAtoms()):
            atom_id = atom.GetSymbol()+ str(atom.GetIdx())
            atom.SetProp("name", atom_id)
        Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

        with open(pkl_file, "wb") as f:
            pickle.dump(mol_sdf, f)
    else:
        if VERBOSE and not DEBUG: print(f"[WARNING] {pkl_file} exists. Will not remake.")

def lookup_compound_id(vault_id, compound_record):
    if 'vault_id' not in compound_record.columns or 'compound_id' not in compound_record.columns:
        print('[ERROR] add vault_id and compound_id columns')
        return None
    match = compound_record.loc[compound_record['vault_id'] == vault_id, 'compound_id']
    compound_id = match.iloc[0] if not match.empty else None
    return compound_id

def unique_ccd(ccd_db, len=5, max_attempts=30):
    '''Generates a unique alphanumeric string of specified length and ensures uniqueness.'''

    chars = string.ascii_uppercase + string.digits

    for _ in range(max_attempts):
        ccd = ''.join(random.choices(chars, k=len))
        if not os.path.exists(f"{ccd_db}/{ccd}.pkl"):
            return ccd
    raise RuntimeError("[ERROR] Could not find a unique CCD ID after max attempts.")
