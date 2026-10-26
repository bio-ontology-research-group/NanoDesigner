import os
import subprocess
import pandas as pd
import time
from Bio import PDB
from shutil import rmtree
import traceback
import sys
from .complex_analysis import clean_extended
sys.path.append("/home/rioszemm/NanoDesigner/dyMEAN") # update
from dyMEAN.data.pdb_utils import Protein, Peptide


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
    This function is to extract the [(chain,position,residue) from a PDB given a specified chain
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
                        res_id = residue_pos_tup[1]
                        res_id = residue.id[1]
                        # res_name = PDB.Polypeptide.three_to_one(residue.get_resname())
                        res_name = three_to_one.get(residue.get_resname())
                        if not residue_pos_tup[2].isalpha():
                            pdb_info.append((chain_id, res_id, res_name))
    return pdb_info


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
                        # res_name = PDB.Polypeptide.three_to_one(residue.get_resname())
                        res_name = three_to_one.get(residue.get_resname())
                        pdb_sequence += res_name
                        if not residue_pos_tup[2].isalpha():
                            pdb_info.append((chain_id, res_id, res_name))

    start_index = pdb_sequence.find(sequence)
    end_index = start_index + len(sequence) - 1

    if start_index == -1:   # -1 if the sequence was not found
        return None

    seq_info =  pdb_info[start_index:end_index + 1]

    return seq_info


def clean_up_files(directory, file_name):
    list_dir = os.listdir(directory)
    for file in list_dir:
        if file_name in file:
            full_path = os.path.join(directory, file)
            try:
                os.remove(full_path)
                # print(f"Deleted file: {full_path}")
            except FileNotFoundError:
                # print(f"File not found, could not delete: {full_path}")
                pass
            except Exception as e:
                print(f"Error deleting file {full_path}: {e}")



def interacting_residues(item, pdb_file, antigen, tmp_dir):

    start_time = time.time() 
    # # Change the current working directory to the output directory
    out_dir = tmp_dir
    current_working_dir = os.getcwd()
    parent_directory = os.path.dirname(current_working_dir)
    dSASA_path = os.path.join(parent_directory, "dr_sasa_n", "build", "dr_sasa")


    command = [dSASA_path, "-m", "4", "-i", pdb_file]

    try:
        # Execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmp_dir)
        stdout, stderr = process.communicate()

        # Check if the command was successful
        if process.returncode != 0:
            print("Command failed!")
            print("Error message:", stderr.decode())
            traceback.print_exc()
            return None  
        else:
            pass
            
    except Exception as e:
        print("An error occurred while executing the command:")
        print(e)
        print(f"Double check full-path command: {dSASA_path}")
        traceback.print_exc()
        return None  



    #Construct the output tsv file from SASA analysis
    folder_path, file_name = os.path.split(pdb_file)
    file_name = os.path.splitext(file_name)[0]
    heavychain= item["heavy_chain"]
    results_tsv_file = os.path.join(out_dir, file_name + f".{heavychain}_vs_{antigen}.by_res.tsv")

    # if the file we need does not exist, we cannot proceed.
    if not os.path.exists(results_tsv_file):
        clean_up_files(out_dir, file_name)
        return None

    # Load the TSV file into a DataFrame
    try:
        df = pd.read_csv(results_tsv_file, sep='\t', index_col=0)
    except Exception as e:
        # if by any reason we create a dataframe i cannot proceed, delete temporal files
        clean_up_files(out_dir, file_name)
        return None

    # Use boolean indexing to extract values different than 0
    non_zero_values = df[df != 0].stack()

    # Check if non_zero_values is empty
    if non_zero_values.empty:
        clean_up_files(out_dir, file_name)
        return None

    # Extract row and column indices of non-zero values
    non_zero_indices = non_zero_values.index.tolist()
    unique_positions = set()

    contacts = []
    contacts_more_info = []
    contacts_complete_info = []
    # Iterate through each non-zero index
    for index in non_zero_indices:
        row, col = index

        # Check if row and col are strings
        if not isinstance(row, str) or not isinstance(col, str):
            continue
        
        row_aa, row_chain, row_pos = row.split('/')
        col_aa, col_chain, col_pos = col.split('/')

        if not row_pos.isnumeric() or not col_pos.isnumeric():
            continue

        if (row_chain, row_pos) in unique_positions or (col_chain, col_pos) in unique_positions:
            continue

        if df.at[row, col] < 1.0:
            continue
    

        row_aa_code = aa_dict.get(row_aa, row_aa)
        col_aa_code = aa_dict.get(col_aa, col_aa)
        
        # print(f"Contact between {row_chain} {row_aa_code}{row_pos} and {col_chain} {col_aa_code}{col_pos}: Value: {df.at[row, col]}")
        contact_tuple = [(row_chain,int(row_pos)), (col_chain,int(col_pos))]
        contacts.append(contact_tuple)

        tup_info_extra = (row_chain, int(row_pos), row_aa_code)
        contacts_more_info.append(tup_info_extra)

        tup_info_extra = (col_chain, int(col_pos), col_aa_code)
        contacts_more_info.append(tup_info_extra)

        contact_tuple_more_info = [(row_chain,int(row_pos),row_aa_code), (col_chain,int(col_pos),col_aa_code)]
        contacts_complete_info.append(contact_tuple_more_info)

        # Add the positions to the set
        unique_positions.add((row_chain, row_pos))
        unique_positions.add((col_chain, col_pos))


    #if the contacts list are empty, stop the analysis, no contacts were found
    if len(contacts) == 0:
        clean_up_files(out_dir, file_name)
        return None

    end_time = time.time() 
    print(f"Interacting aminoacids computation took: {end_time - start_time:.2f} seconds")

    nanobody_aa_contact = []
    for element in contacts_more_info:
        if element[0] == heavychain:  
            nanobody_aa_contact.append(element)

    antigen_aa_contact = []
    for element in contacts_more_info:
        if element[0] == antigen:  
            antigen_aa_contact.append(element)
    

    if not "cdrh3_seq_mod" in item:
        cdrh3_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh3_seq"])
        cdrh2_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh2_seq"])
        cdrh1_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh1_seq"])
    else:
        cdrh3_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh3_seq_mod"])
        cdrh2_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh2_seq_mod"])
        cdrh1_aminoacids = extract_seq_info_from_pdb(pdb_file, heavychain, item["cdrh1_seq_mod"])  


    # Check if any of the CDR lists are None and set them to empty lists if they are
    if cdrh3_aminoacids is None:
        print("cdrh3_aminoacids is None, setting to empty list")
        cdrh3_aminoacids = []
    if cdrh2_aminoacids is None:
        print("cdrh2_aminoacids is None, setting to empty list")
        cdrh2_aminoacids = []
    if cdrh1_aminoacids is None:
        print("cdrh1_aminoacids is None, setting to empty list")
        cdrh1_aminoacids = []

    cdr3_matches = 0
    cdr3_matching_res = []
    for element in cdrh3_aminoacids:
        if element in nanobody_aa_contact:
            cdr3_matches += 1
            cdr3_matching_res.append(element)
            
    cdr2_matches = 0
    cdr2_matching_res = []
    for element in cdrh2_aminoacids:
        if element in nanobody_aa_contact:
            cdr2_matches += 1
            cdr2_matching_res.append(element)

    cdr1_matches = 0
    cdr1_matching_res = []
    for element in cdrh1_aminoacids:
        if element in nanobody_aa_contact:
            cdr1_matches += 1
            cdr1_matching_res.append(element)



    cdr3_involvement = (cdr3_matches/len(cdrh3_aminoacids))*100 if len(cdrh3_aminoacids) != 0.0 else 0.0
    cdr2_involvement = (cdr2_matches/len(cdrh2_aminoacids))*100 if len(cdrh2_aminoacids) != 0.0 else 0.0
    cdr1_involvement = (cdr1_matches/len(cdrh1_aminoacids))*100 if len(cdrh1_aminoacids) != 0.0 else 0.0


    # Print the average involvement of each CDR
    # print("Average involvement of CDR H3:", cdr3_involvement)
    # print("Average involvement of CDR H2:", cdr2_involvement)
    # print("Average involvement of CDR H1:", cdr1_involvement)
    
    # item["cdr3_avg"] = cdr3_involvement
    # item["cdr2_avg"] = cdr2_involvement
    # item["cdr1_avg"] = cdr1_involvement

    # Define epitope as antigen residues in touch with the cdrs
    antigen_aa_contact = [(tuple_item[0], tuple_item[1]) for tuple_item in antigen_aa_contact]

    cdrs_aa = cdrh3_aminoacids[:]
    cdrs_aa.extend(cdrh2_aminoacids)
    cdrs_aa.extend(cdrh1_aminoacids)


    # Extract epitope
    contact_dict = {tuple(contact[1]): tuple(contact[0]) for contact in contacts_complete_info}
    first_key = next(iter(contact_dict))
    
    if first_key[0] != heavychain:
        contact_dict = swapped_dict = {v: k for k, v in contact_dict.items()}

    keys_list = list(contact_dict.keys())

    epitope = []
    for tup in cdrs_aa:
        if tup in keys_list:
            epi_aa = contact_dict[tup]
            epitope.append(epi_aa)


    clean_up_files(out_dir, file_name)

    return epitope, cdr3_matching_res, cdr2_matching_res, cdr1_matching_res


def interacting_residues_extended(item, pdb_file, antigen, tmp_dir, cdr_type:None):

    # by defaul the CDRHs
    if cdr_type == None:
        cdr_type = "H"
    # else L 

    start_time = time.time() 
    # # Change the current working directory to the output directory
    out_dir = tmp_dir
    current_working_dir = os.getcwd()
    parent_directory = os.path.dirname(current_working_dir)
    dSASA_path = os.path.join(parent_directory, "dr_sasa_n", "build", "dr_sasa")
    command = [dSASA_path, "-m", "4", "-i", pdb_file]

    try:
        # Execute the command
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tmp_dir)
        stdout, stderr = process.communicate()

        # Check if the command was successful
        if process.returncode != 0:
            print("Command failed!")
            print("Error message:", stderr.decode())
            traceback.print_exc()
            return None  
        else:
            print("Command executed successfully.")
            
    except Exception as e:
        print("An error occurred while executing the command:")
        print(e)
        print(f"Double check full-path command: {dSASA_path}")
        traceback.print_exc()
        return None


    #Construct the output tsv file from SASA analysis
    folder_path, file_name = os.path.split(pdb_file)
    file_name = os.path.splitext(file_name)[0]
    heavychain= item["heavy_chain"]
    lightchain= item.get("light_chain", None)

    if cdr_type == "H":
        results_tsv_file = os.path.join(out_dir, file_name + f".{heavychain}_vs_{antigen}.by_res.tsv")
    else:
        results_tsv_file = os.path.join(out_dir, file_name + f".{lightchain}_vs_{antigen}.by_res.tsv")

    # if the file we need does not exist, we cannot proceed.
    if not os.path.exists(results_tsv_file):
        clean_up_files(out_dir, file_name)
        return None

    # Load the TSV file into a DataFrame
    try:
        df = pd.read_csv(results_tsv_file, sep='\t', index_col=0)
    except Exception as e:
        # if by any reason we create a dataframe i cannot proceed, delete temporal files
        clean_up_files(out_dir, file_name)
        return None

    # Use boolean indexing to extract values different than 0
    non_zero_values = df[df != 0].stack()

    # Check if non_zero_values is empty
    if non_zero_values.empty:
        clean_up_files(out_dir, file_name)
        return None

    # Extract row and column indices of non-zero values
    non_zero_indices = non_zero_values.index.tolist()
    unique_positions = set()

    contacts = []
    contacts_more_info = []
    contacts_complete_info = []
    # Iterate through each non-zero index
    for index in non_zero_indices:
        row, col = index

        # Check if row and col are strings
        if not isinstance(row, str) or not isinstance(col, str):
            continue
        
        row_aa, row_chain, row_pos = row.split('/')
        col_aa, col_chain, col_pos = col.split('/')

        if not row_pos.isnumeric() or not col_pos.isnumeric():
            continue

        if (row_chain, row_pos) in unique_positions or (col_chain, col_pos) in unique_positions:
            continue

        if df.at[row, col] < 1.0:
            continue
    

        row_aa_code = aa_dict.get(row_aa, row_aa)
        col_aa_code = aa_dict.get(col_aa, col_aa)
        
        # print(f"Contact between {row_chain} {row_aa_code}{row_pos} and {col_chain} {col_aa_code}{col_pos}: Value: {df.at[row, col]}")
        contact_tuple = [(row_chain,int(row_pos)), (col_chain,int(col_pos))]
        contacts.append(contact_tuple)

        tup_info_extra = (row_chain, int(row_pos), row_aa_code)
        contacts_more_info.append(tup_info_extra)

        tup_info_extra = (col_chain, int(col_pos), col_aa_code)
        contacts_more_info.append(tup_info_extra)

        contact_tuple_more_info = [(row_chain,int(row_pos),row_aa_code), (col_chain,int(col_pos),col_aa_code)]
        contacts_complete_info.append(contact_tuple_more_info)

        # Add the positions to the set
        unique_positions.add((row_chain, row_pos))
        unique_positions.add((col_chain, col_pos))


    #if the contacts list are empty, stop the analysis, no contacts were found
    if len(contacts) == 0:
        clean_up_files(out_dir, file_name)
        return None

    end_time = time.time() 
    print(f"Interacting aminoacids computation took: {end_time - start_time:.2f} seconds")

    chain_key = {"H": "heavy_chain", "L": "light_chain"}.get(cdr_type)
    chain = item.get(chain_key)

    nanobody_aa_contact = []
    for element in contacts_more_info:
        if element[0] == chain:  
            nanobody_aa_contact.append(element)

    antigen_aa_contact = []
    for element in contacts_more_info:
        if element[0] == antigen:  
            antigen_aa_contact.append(element)

    
    # Define keys dynamically based on `cdr_type`
    cdr_keys = {
        "cdr3_seq": f"cdr{cdr_type.lower()}3_seq",
        "cdr2_seq": f"cdr{cdr_type.lower()}2_seq",
        "cdr1_seq": f"cdr{cdr_type.lower()}1_seq",
        "cdr3_seq_mod": f"cdr{cdr_type.lower()}3_seq_mod",
        "cdr2_seq_mod": f"cdr{cdr_type.lower()}2_seq_mod",
        "cdr1_seq_mod": f"cdr{cdr_type.lower()}1_seq_mod"
    }

    # Extract CDR sequences, with modified sequences if available
    cdr3_seq = item.get(cdr_keys["cdr3_seq_mod"], item[cdr_keys["cdr3_seq"]])
    cdr2_seq = item.get(cdr_keys["cdr2_seq_mod"], item[cdr_keys["cdr2_seq"]])
    cdr1_seq = item.get(cdr_keys["cdr1_seq_mod"], item[cdr_keys["cdr1_seq"]])

    cdr3_aminoacids = extract_seq_info_from_pdb(pdb_file, chain, cdr3_seq)
    cdr2_aminoacids = extract_seq_info_from_pdb(pdb_file, chain, cdr2_seq)
    cdr1_aminoacids = extract_seq_info_from_pdb(pdb_file, chain, cdr1_seq)


    # Ensure amino acid lists are not None
    cdr3_aminoacids = cdr3_aminoacids if cdr3_aminoacids is not None else []
    cdr2_aminoacids = cdr2_aminoacids if cdr2_aminoacids is not None else []
    cdr1_aminoacids = cdr1_aminoacids if cdr1_aminoacids is not None else []

    # # Calculate matches and involvement for each CDR
    # def calculate_matches_and_involvement(cdr_aminoacids):
    #     matches = 0
    #     matching_res = []
    #     for element in cdr_aminoacids:
    #         if element in nanobody_aa_contact:
    #             matches += 1
    #             matching_res.append(element)
    #     involvement = (matches / len(cdr_aminoacids)) * 100 if cdr_aminoacids else 0.0
    #     return matches, matching_res, involvement

    # cdr3_matches, cdr3_matching_res, cdr3_involvement = calculate_matches_and_involvement(cdr3_aminoacids)
    # cdr2_matches, cdr2_matching_res, cdr2_involvement = calculate_matches_and_involvement(cdr2_aminoacids)
    # cdr1_matches, cdr1_matching_res, cdr1_involvement = calculate_matches_and_involvement(cdr1_aminoacids)


    cdrs_aa = cdr3_aminoacids + cdr2_aminoacids + cdr1_aminoacids 
    contact_dict = {tuple(contact[1]): tuple(contact[0]) for contact in contacts_complete_info}
    first_key = next(iter(contact_dict))
    
    if first_key[0] != chain:
        contact_dict = swapped_dict = {v: k for k, v in contact_dict.items()}

    keys_list = list(contact_dict.keys())


    epitope = []
    for tup in cdrs_aa:
        if tup in keys_list:
            epi_aa = contact_dict[tup]
            epitope.append(epi_aa)

    return epitope 

# 1
def binding_residues_analysis(shared_args,hdock_model):

    print("conducting binding_residues_analysis function")

    start_time = time.time()
    ag_pdb, ab_pdb, binding_rsite, item, pdb_dir_hdock, tmp_dir_for_interacting_aa = shared_args

    antigen_chains = item["antigen_chains"]

    cdrh3_seq = item["cdrh3_seq_mod"]
    cdrh2_seq = item["cdrh2_seq_mod"]
    cdrh1_seq = item["cdrh1_seq_mod"] 

    # sanity check, make all the models to have same numbering as the ground truth pdbs
    H, L, A = item['heavy_chain'], item['light_chain'], item['antigen_chains']
    mod_revised, ref_revised = False, False
    chains_list = [H] + A
    if L != "":
        chains_list = [H] + A + [L]

    try:
        mod_prot = Protein.from_pdb(hdock_model,chains_list)
    except Exception as e:
        print(f'parse {hdock_model} failed for {e}')

    try:
        ref_prot = Protein.from_pdb(ab_pdb,chains_list)
    except Exception as e:
        print(f'parse {hdock_model} failed for {e}')

    try:
        # Iterate through chains specified in the list
        for chain_name in chains_list:
            # Retrieve Peptide objects by chain name
            ref_chain = ref_prot.peptides.get(chain_name)
            mod_chain = mod_prot.peptides.get(chain_name)

            # Check if chains exist in both structures
            if mod_chain is None:
                print(f"Chain {chain_name} not found in one of the complexes. Skipping.")
                continue

            if ref_chain is None:
                pass


            # Align chains by position number if their lengths differ
            if len(mod_chain.residues) != len(ref_chain.residues):
                print(f"{mod_chain} chain {chain_name} length not consistent: {len(mod_chain.residues)} != {len(ref_chain.residues)}. Trying to align by position number.")

                # Prepare to store revised residues
                mod_residues, ref_residues = [], []
                pos_map = {'-'.join(str(a) for a in res.get_id()): res for res in ref_chain.residues}

                for res in mod_chain.residues:
                    res_id = '-'.join(str(a) for a in res.get_id())
                    if res_id in pos_map:
                        mod_residues.append(res)
                    else:
                        # Handle the case where the residue position does not match
                        print(f"Warning: Residue {res_id} in mod_chain not found in ref_chain. Skipping.")

                ref_residues = [pos_map[res_id] for res_id in [ '-'.join(str(a) for a in res.get_id()) for res in mod_residues ] if res_id in pos_map]

                # Update peptide chains within the Protein objects 
                mod_prot.peptides[chain_name] = Peptide(chain_name, mod_residues)
                ref_prot.peptides[chain_name] = Peptide(chain_name, ref_residues)

                mod_revised, ref_revised = True, False  # Only mod_prot is revised
                print(f"{mod_prot.peptides[chain_name]} chain {chain_name} length after aligned: {len(mod_prot.peptides[chain_name].residues)} == {len(ref_prot.peptides[chain_name].residues)}")

    except Exception as e:
        print(f"An exception was raised: {e}")

    if mod_revised:
        print("Entered the mod_revised mode")
        mod_prot.to_pdb(hdock_model)


    try:
        clean_extended(ag_pdb, ab_pdb, hdock_model, hdock_model)
    except Exception as e:
        print(f"model {hdock_model}, gave this error: {e}")
        print(traceback.format_exc())

    # The next step is to check the involvement of both, ground truth epitope and paratope
    gt_epitope = binding_rsite  # epitope, which comes from the ground truth

    # Get the base name of the file
    filename = os.path.basename(hdock_model)

    # Split the filename into name and extension
    model_name, extension = os.path.splitext(filename)

    epitope_model = []
    cdr1_interactions_to_ag = []
    cdr2_interactions_to_ag = []
    cdr3_interactions_to_ag = []

    print(f"Starting time before binding interface analysis")
    start_time1 = time.time()
    for antigen in antigen_chains:
        print(f"Processing {antigen}")
        # the program doesnt work with PDBs with one more than one chain in the PDB; create a temp pdb with the current pdb
        if len(antigen_chains) > 1:
            tmp_pdb = os.path.join(tmp_dir_for_interacting_aa, f"{model_name}_chain{antigen}.pdb")
            chains_to_reconstruct = []
            chains_to_reconstruct.extend(H)
            chains_to_reconstruct.extend(antigen)
            # program crashes with more than 2 chains, as we only care about cdrh involvement we will skip reconstruction of light chain in this temporal file
            # if L:
            #     chains_to_reconstruct.extend(L)
            try:
                protein = Protein.from_pdb(hdock_model, chains_to_reconstruct)
                protein.to_pdb(tmp_pdb)
            except Exception as e:
                print(f"Failed to process PDB file '{hdock_model}': {e}")
                continue

            result = interacting_residues(item, tmp_pdb, antigen, tmp_dir_for_interacting_aa)
            if result is not None:
                # epitope_result, cdr3_matching_res, cdr2_matching_res, cdr1_matching_res = interacting_residues(item, tmp_pdb, antigen,tmp_dir)
                epitope_result, cdr3_matching_res, cdr2_matching_res, cdr1_matching_res = result
            else:
                continue
        else:
            result = interacting_residues(item, hdock_model, antigen, tmp_dir_for_interacting_aa)
            if result is not None:
                # epitope_result, cdr3_matching_res, cdr2_matching_res, cdr1_matching_res = interacting_residues(item, tmp_pdb, antigen,tmp_dir)
                epitope_result, cdr3_matching_res, cdr2_matching_res, cdr1_matching_res = result
            else:
                continue
        
        end_time1 = time.time()
        print(f"Finish binding residues analysis on {antigen} : {end_time1 - start_time1:.2f} seconds")

        epitope_model.extend(epitope_result) 
        cdr3_interactions_to_ag.extend(cdr3_matching_res)
        cdr2_interactions_to_ag.extend(cdr2_matching_res)
        cdr1_interactions_to_ag.extend(cdr1_matching_res)

    cdr1_interactions_to_ag = [list(t) for t in set(tuple(i) for i in cdr1_interactions_to_ag)]
    cdr2_interactions_to_ag = [list(t) for t in set(tuple(i) for i in cdr2_interactions_to_ag)]
    cdr3_interactions_to_ag = [list(t) for t in set(tuple(i) for i in cdr3_interactions_to_ag)]


    cdr3_involvement = (len(cdr3_interactions_to_ag)/len(cdrh3_seq))*100 if len(cdrh3_seq) != 0.0 else 0.0
    cdr2_involvement = (len(cdr2_interactions_to_ag)/len(cdrh2_seq))*100 if len(cdrh2_seq) != 0.0 else 0.0
    cdr1_involvement = (len(cdr1_interactions_to_ag)/len(cdrh1_seq))*100 if len(cdrh1_seq) != 0.0 else 0.0
    total_avg_cdr_involvement = float((cdr3_involvement + cdr2_involvement + cdr1_involvement) / 3)

    epitope_model = [tuple(element) for element in epitope_model]
    epitope_model = [tuple(element[:2]) for element in epitope_model]
    epitope_model = list(set(epitope_model))

    # how much of the gt_epitope was recovered?
    # Convert lists of tuples to sets
    epitope_model_set = set(epitope_model)
    gt_epitope_set = set(gt_epitope)

    # Find the common tuples (order may not be preserved)
    common_amino_acids_set = epitope_model_set.intersection(gt_epitope_set)
    common_amino_acids_list = list(common_amino_acids_set)

    # Calculate recall
    epitope_recall = len(common_amino_acids_list) / len(gt_epitope)

    end_time = time.time() 
    print(f"Process_hdock_model took: {end_time - start_time:.2f} seconds")

    tup_result = (hdock_model, cdr1_involvement, cdr2_involvement, cdr3_involvement, total_avg_cdr_involvement, epitope_recall)
    keys = ['hdock_model', 'cdr1_avg', 'cdr2_avg', 'cdr3_avg', 'total_avg_cdr_involvement', 'epitope_recall']
    result_dict = {keys[i]: tup_result[i] for i in range(len(keys))}

    return result_dict

