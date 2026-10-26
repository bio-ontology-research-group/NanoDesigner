import sys
import numpy as np
sys.path.append("/home/rioszemm/NanoDesigner/dyMEAN") # update
from dyMEAN.data.pdb_utils import VOCAB, Protein
from dyMEAN.configs import IMGT 
from dyMEAN.utils.network import url_get
from Bio import PDB
from scipy.spatial import cKDTree


parser = PDB.PDBParser(QUIET=True)
downloader = PDB.PDBList()



def fetch_from_sabdab(identifier, numbering_scheme, save_path, tries=5):
    print(f"fetching {identifier} from sabdab")
    
    identifier = identifier.lower()
    url = f'https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/pdb/{identifier}/?scheme={numbering_scheme}'

    text = url_get(url, tries)
    if text is None:
        return None

    with open(save_path, 'w') as file:
        file.write(text.text)

    data = {
        'pdb': save_path
    }
    return data

def fetch_from_pdb(identifier, save_path, tries=5):
    # example identifier: 1FBI

    identifier = identifier.upper()
    url = 'https://data.rcsb.org/rest/v1/core/entry/' + identifier

    res = url_get(url, tries)
    if res is None:
        return None

    url = f'https://files.rcsb.org/download/{identifier}.pdb'


    text = url_get(url, tries)
    if text is None:
        return None

    with open(save_path, 'w') as file:
        file.write(text.text)

    data = {
        'pdb': save_path
    }


def extract_antibody_info(antibody, heavy_ch, light_ch, numbering):

    _scheme = IMGT
    if numbering == 'imgt':
        _scheme = IMGT

    # get cdr/frame denotes
    h_type_mapping, l_type_mapping = {}, {}  # - for non-Fv region, 0 for framework, 1/2/3 for cdr1/2/3

    for lo, hi in [_scheme.HFR1, _scheme.HFR2, _scheme.HFR3, _scheme.HFR4]:
        for i in range(lo, hi + 1):
            h_type_mapping[i] = '0'
    for cdr, (lo, hi) in zip(['1', '2', '3'], [_scheme.H1, _scheme.H2, _scheme.H3]):
        for i in range(lo, hi + 1):
            h_type_mapping[i] = cdr
    h_conserved = _scheme.Hconserve

    for lo, hi in [_scheme.LFR1, _scheme.LFR2, _scheme.LFR3, _scheme.LFR4]:
        for i in range(lo, hi + 1):
            l_type_mapping[i] = '0'
    for cdr, (lo, hi) in zip(['1', '2', '3'], [_scheme.L1, _scheme.L2, _scheme.L3]):
        for i in range(lo, hi + 1):
            l_type_mapping[i] = cdr
    l_conserved = _scheme.Lconserve

    # get variable domain and cdr positions
    selected_peptides, cdr_pos = {}, {}

    chain_names = [heavy_ch]
    if light_ch:
        chain_names.append(light_ch)

    for c, chain_name in zip(['H', 'L'], chain_names):
        chain = antibody.get_chain(chain_name)
        # print("chain_names",chain_names)
        if chain is None:
            continue  # Skip processing if the chain is None
        # Note: possbly two chains are different segments of a same chain
        assert chain is not None, f'Chain {chain_name} not found in the antibody'
        type_mapping = h_type_mapping if c == 'H' else l_type_mapping
        conserved = h_conserved if c == 'H' else l_conserved
        res_type = ''
        for i in range(len(chain)):
            residue = chain.get_residue(i)
            residue_number = residue.get_id()[0]

            if residue_number in type_mapping:
                res_type += type_mapping[residue_number]
                if residue_number in conserved:
                    hit, symbol = False, residue.get_symbol()
                    for conserved_residue in conserved[residue_number]:
                        if symbol == VOCAB.abrv_to_symbol(conserved_residue):
                            hit = True
                            break
                        print(f"Actual residue type at position {residue_number} in chain {chain_name}: {symbol}")
                    assert hit, f'Not {conserved[residue_number]} at {residue_number} in chain {chain_name}'
            else:
                res_type += '-'
        if '0' not in res_type:
            print(heavy_ch, light_ch, antibody.pdb_id, res_type)
        start, end = res_type.index('0'), res_type.rindex('0')
        for cdr in ['1', '2', '3']:
            cdr_start, cdr_end = res_type.find(cdr), res_type.rfind(cdr)
            assert cdr_start != -1, f'cdr {c}{cdr} not found, residue type: {res_type}'
            start, end = min(start, cdr_start), max(end, cdr_end)
            cdr_pos[f'CDR-{c}{cdr}'] = (cdr_start, cdr_end)
        for cdr in ['1', '2', '3']:
            cdr = f'CDR-{c}{cdr}'
            cdr_start, cdr_end = cdr_pos[cdr]
            cdr_pos[cdr] = (cdr_start - start, cdr_end - start)
        chain = chain.get_span(start, end + 1)  # the length may exceed 130 because of inserted amino acids
        chain.set_id(chain_name)
        selected_peptides[chain_name] = chain

    return cdr_pos

def get_cdr_pos(cdr_pos, cdr):  # H/L + 1/2/3, return [begin, end] position
    # cdr_pos is a dictionary = {'CDR-H1': (24, 31), 'CDR-H2': (49, 56), 'CDR-H3': (95, 114)}
    cdr = f'CDR-{cdr}'.upper()
    if cdr in cdr_pos:
        return cdr_pos[cdr]
    else:
        return None


def extract_seq_info_from_pdb(pdb_file, target_chain,sequence):
    """
    This function is to extract the [(chain,position,residue) from a PDB given a specified chain
    """
    structure = PDB.PDBParser().get_structure("protein", pdb_file)
    pdb_info = []
    pdb_sequence = ""
    for model in structure:
        for chain in model:
            if chain.id == target_chain:
                for residue in chain:
                    # residue <Residue VAL het=  resseq=124 icode= >
                    if PDB.is_aa(residue):
                        chain_id = chain.id
                        residue_pos_tup = residue.id
                        #residue_pos_tup: 
                        # (' ', 111, ' ')
                        # (' ', 111, 'A')
                        # (' ', 111, 'B')
                        res_id = residue_pos_tup[1]
                        res_name = PDB.Polypeptide.three_to_one(residue.get_resname())
                        pdb_sequence += res_name
                        if not residue_pos_tup[2].isalpha():
                            pdb_info.append((chain_id, res_id, res_name))


    start_index = pdb_sequence.find(sequence)
    end_index = start_index + len(sequence) - 1

    if start_index == -1:   # -1 if the sequence was not found
        print("Specified sequence was not found in the PDB")
        return None

    seq_info =  pdb_info[start_index:end_index + 1]

    return seq_info



aa_dict = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 
    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 
    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N', 
    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S', 
    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_pdb_info(pdb_file, target_chain):
    """
    Extract [(chain, position, residue)] from a PDB given a specified chain.
    Handles insertion codes (like 100A, 100B) properly.
    """
    structure = PDB.PDBParser().get_structure("protein", pdb_file)
    pdb_info = []
    for model in structure:
        for chain in model:
            if chain.id == target_chain:
                for residue in chain:
                    if PDB.is_aa(residue):
                        chain_id = chain.id
                        residue_pos_tup = residue.id
                        # Handle residue numbering with insertion codes
                        res_id = residue_pos_tup[1]
                        insertion_code = residue_pos_tup[2].strip()  # Get insertion code if exists
                        if insertion_code:
                            res_id = f"{res_id}{insertion_code}"  # Combine number and insertion code
                        
                        res_name = three_to_one.get(residue.get_resname(), 'X')  # Default to 'X' if unknown
                        pdb_info.append((chain_id, res_id, res_name))
    return pdb_info


    

def clean_extended(origin_antigen_pdb, origin_antibody_pdb, template_pdb, out_pdb):
    origin_antigen_cplx = Protein.from_pdb(origin_antigen_pdb)
    origin_antibody_cplx = Protein.from_pdb(origin_antibody_pdb)
    template_cplx = Protein.from_pdb(template_pdb)
    peptides = {}

    ori_antigen_chain_to_id, ori_antibody_chain_to_id = {}, {}
    id_to_temp_chain = {}
    
    for chain_name, chain in origin_antigen_cplx:
        ori_antigen_chain_to_id[chain_name] = f'{chain.get_seq()[:5]}'
    for chain_name, chain in origin_antibody_cplx:
        ori_antibody_chain_to_id[chain_name] = f'{chain.get_seq()[:5]}'
    for chain_name, chain in template_cplx:
        id_to_temp_chain[f'{chain.get_seq()[:5]}'] = chain_name
    
    for chain_name in origin_antigen_cplx.get_chain_names():
        ori_chain = origin_antigen_cplx.get_chain(chain_name)
        temp_chain = template_cplx.get_chain(id_to_temp_chain[ori_antigen_chain_to_id[chain_name]])
        for i, residue in enumerate(ori_chain):
            if i < len(temp_chain):
                # renumber
                temp_chain.residues[i].id = residue.id
                # delete Hs
                for atom in temp_chain.residues[i].coordinate:
                    if atom[0] == 'H':
                        del temp_chain.residues[i].coordinate[atom]
            else:
                print(f'{origin_antigen_cplx.get_id()}, chain {chain_name} lost residues {len(ori_chain)} > {len(temp_chain)}')
                break
        temp_chain.set_id(chain_name)
        peptides[chain_name] = temp_chain
    
    for chain_name in origin_antibody_cplx.get_chain_names():
        ori_chain = origin_antibody_cplx.get_chain(chain_name)
        temp_chain = template_cplx.get_chain(id_to_temp_chain[ori_antibody_chain_to_id[chain_name]])
        for i, residue in enumerate(ori_chain):
            if i < len(temp_chain):
                # renumber
                temp_chain.residues[i].id = residue.id
                # delete Hs
                for atom in temp_chain.residues[i].coordinate:
                    if atom[0] == 'H':
                        del temp_chain.residues[i].coordinate[atom]
            else:
                print(f'{origin_antibody_cplx.get_id()}, chain {chain_name} lost residues {len(ori_chain)} > {len(temp_chain)}')
                break
        temp_chain.set_id(chain_name)
        peptides[chain_name] = temp_chain
    
    renumber_cplx = Protein(template_cplx.get_id(), peptides)
    renumber_cplx.to_pdb(out_pdb)



# Atomic radii for various atom types. 
# You can comment out the ones you don't care about or add new ones
atom_radii = {
#    "H": 1.20,  # Who cares about hydrogen??
    "C": 1.70, 
    "N": 1.55, 
    "O": 1.52,
    "S": 1.80,
    "F": 1.47, 
    "P": 1.80, 
    "CL": 1.75, 
    "MG": 1.73,
}

def count_clashes(structure, clash_cutoff=0.63):

    """
    SOURCE: https://www.blopig.com/blog/2023/05/checking-your-pdb-file-for-clashing-atoms/
    with additional information about the clashes
    """

    # Set what we count as a clash for each pair of atoms
    clash_cutoffs = {i + "_" + j: (clash_cutoff * (atom_radii[i] + atom_radii[j])) for i in atom_radii for j in atom_radii}

    # Extract atoms for which we have radii
    atoms = [x for x in structure.get_atoms() if x.element in atom_radii]
    coords = np.array([a.coord for a in atoms], dtype="d")

    # Build a KDTree using scipy
    kdt = cKDTree(coords)

    # Initialize a list to hold clashes and clash_details
    clashes = []
    clash_details = []

    # Iterate through all atoms
    for atom_1 in atoms:
        # Find atoms that could be clashing
        kdt_search = kdt.query_ball_point(atom_1.coord, max(clash_cutoffs.values()))

        # Get index and distance of potential clashes
        potential_clash = [(ix, np.linalg.norm(coords[ix] - atom_1.coord)) for ix in kdt_search]

        for ix, atom_distance in potential_clash:
            atom_2 = atoms[ix]

            # Exclude clashes from atoms in the same residue, peptide bonds, and disulphide bridges (similar to your original code)
            if (
                atom_1.parent.id == atom_2.parent.id
                or (atom_2.name == "C" and atom_1.name == "N")
                or (atom_2.name == "N" and atom_1.name == "C")
                or (atom_2.name == "SG" and atom_1.name == "SG" and atom_distance > 1.88)
            ):
                continue

            clash_type = f"{atom_1.element}-{atom_2.element} clash"
            if atom_distance < clash_cutoffs[atom_2.element + "_" + atom_1.element]:
                clashes.append((atom_1, atom_2))
                clash_details.append({
                    'clash_type': clash_type,
                    'atoms_involved': (atom_1, atom_2),
                    'distance': atom_distance,
                    'residues_involved': (atom_1.parent.id, atom_2.parent.id)
                })

    return len(clashes) // 2, list(clash_details)
