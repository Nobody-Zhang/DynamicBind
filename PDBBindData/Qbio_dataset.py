import re
import requests
import pandas as pd
import os
import threading

import shutil

from gemmi import cif
import shutil
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from rdkit import Chem
import concurrent.futures

from tqdm import tqdm
import numpy as np

df_result = pd.DataFrame()
lock = threading.Lock()
lock2 = threading.Lock()

# from openbabel.pybel import (readfile,Outputfile) 

def MolFormatConversion(input_file, output_file):
    shutil.copy2(input_file, output_file)
    # mol = Chem.MolFromPDBFile(input_file)
    # Chem.MolTomol2File(mol, output_file)

def download_cif(pdb_code, cif_path):
    # Download the cif file
    url = f"https://files.rcsb.org/download/{pdb_code}.cif"
    response = requests.get(url)
    with open(cif_path, 'wb') as f:
        f.write(response.content)

def find_assembly_index(pdb_code, pdb_ori):
    # Read the cif file
    cif_path = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/.tmp/{pdb_code}.cif'
    if not os.path.exists(cif_path):
        download_cif(pdb_code, cif_path)
    doc = cif.read_file(cif_path)
    block = doc.sole_block()

    # Find assembly Chains
    list_assembly_all = list(block.find_values('_pdbx_struct_assembly_gen.asym_id_list'))
    list_assembly_all = [_.split(',') for _ in list_assembly_all]
    checked = 0
    flag = False
    if len(list_assembly_all) > 1:
        # Means that there are multiple assemblies,
        # Need information from the pdbbind website
        sequence_url = f"https://www.rcsb.org/annotations/{pdb_code}"
        resp_seq = requests.get(sequence_url)
        # The annotations website guarantees that there is only protein information
        pattern_auth = re.compile(r'\[auth')
        matches_auth = pattern_auth.finditer(resp_seq.text)

        positions = [match.start() for match in matches_auth]
        dic_name = {}
        if positions != []:
            # This means that there are changes of the chain name
            for _ in positions:
                pdb_name = resp_seq.text[_+6]
                cif_name = resp_seq.text[_-2]
                if pdb_name not in dic_name and cif_name >= 'A' and cif_name <= 'Z':
                    dic_name[pdb_name] = cif_name
        with open(f"/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_ori}/{pdb_ori}_protein_processed.pdb", "r") as f:
            # Read the first line which contains chains selected
            lines = f.readlines()
            # Choose the 4th one, if multiple words, then leave the 1st one
            chains = [line.split()[4][0] for line in lines if line.startswith('ATOM')]
            chains = list(set(chains))
            

            # Now we know the chains contained in the pdb file
            for idx, chain in enumerate(chains):
                if chain in dic_name:
                    chains[idx] = dic_name[chain]
            for idx, chain_list in enumerate(list_assembly_all):
                flag = True
                for chain in chains:
                    if chain not in chain_list:
                        flag = False
                        break
                if flag:
                    # Find the matched chain
                    checked = idx
                    break
    
    # if not flag and len(list_assembly_all) > 1:
    #     print(f"Warning: {pdb_code} has no matched assembly chains, use the first one")
    return checked

cnt = 0
multi_cnt = 0

def handle(pdb, lig):
    global cnt
    global df_result
    global multi_cnt
    global lock
    global lock2
    global multi_bs

    pdb_ori = pdb
    if pdb in removed:
        pdb = removed[pdb]
    dest_dir = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_Qbio/{pdb}'
    filtered_df = df[(df['PDB ID'] == pdb) & (df['Ligand ID'] == lig)]

    os.makedirs(dest_dir, exist_ok=True)
    
    assembly_index = find_assembly_index(pdb, pdb_ori) + 1
    filtered_df = filtered_df[filtered_df['Assembly ID'] == f'{pdb}_{assembly_index}']

    dest_pdb_file = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_Qbio/{pdb}/{pdb}_protein.pdb'

    check_bs = []

    for index, row in filtered_df.iterrows():
        ligand_detail = row['Ligand Detail']
        file = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/lig_pdb/{ligand_detail}.pdb'
        # print(os.path.isfile(file))
        dest_file = os.path.join(dest_dir, f'{ligand_detail}.pdb')
        MolFormatConversion(file, dest_file)

        if row['Site'] not in check_bs:
            check_bs.append(row['Site'])
        else:
            # Same biding site but multiple ligand forms
            with lock2:
                multi_bs.append(pdb)
        
        with lock:
            df_result = df_result.append({'PDB': pdb, 'PDB_file': dest_pdb_file, 'Ligand': lig, 'Ligand_file':dest_file, 'BS': row['Site'], 'Gen': 1}, ignore_index=True)
        # with open(f'{dest_dir}/{ligand_detail}.pdb', 'w') as f:
        #     f.write(resp_lig.text)
    
    if len(filtered_df) > 1: 
        with lock2:
            multi_cnt += 1

    if len(filtered_df) == 0: 
        # probably due to ligand name mismatched
        with lock:
            src_file = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_ori}/{pdb_ori}_ligand.mol2'
            dst_file = os.path.join(dest_dir, f'{lig}.mol2')
            shutil.copy2(src_file, dst_file)
            df_result = df_result.append({'PDB': pdb, 'PDB_file': dest_pdb_file, 'Ligand': lig, 'Ligand_file':dst_file, 'BS': None, 'Gen': 0}, ignore_index=True)
    
    source_pdb_file = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_ori}/{pdb_ori}_protein_processed.pdb'
    shutil.copy2(source_pdb_file, dest_pdb_file)
    with lock2:
        cnt += 1
    # print(f"Finished {cnt / len(pdb_codes) * 100:.2f}% MULTI: {multi_cnt / cnt * 100:.2f}%")
    # print(df_result)

from Bio.PDB import PDBParser

def calculate_distance(ligand_path, protein_path):
    # Load ligand and protein structures
    ligand = Chem.MolFromPDBFile(ligand_path)

    # Parse the PDB file
    parser = PDBParser()
    structure = parser.get_structure('protein', protein_path)

    # Get the coordinates of each amino acid
    coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())
    
    # Get atom coordinates for ligand and protein
    ligand_coords = np.array(ligand.GetConformer(0).GetPositions())

    protein_coords = np.array(coordinates)

    # Calculate distances between each ligand atom and protein atom
    distances = np.linalg.norm(ligand_coords[:, np.newaxis] - protein_coords, axis=2)

    # Check if any distance is less than 5A
    if np.any(distances < 5):
        return False
    else:
        return True


if __name__ == '__main__':
    # Thread pool
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=(os.cpu_count()))
    # dictionary to store the mapping information from pdb_code to ligand_name, extracted from INDEX_general_PL_data.2020
    pdbbind_index = {}
    # Parse PDBBind index
    with open('/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/INDEX_general_PL_data.2020', 'r') as f:
        # Skip the header lines
        for _ in range(6):
            next(f)
        
        # Parse the rest of the lines
        for line in f:
            # Split the line into columns
            columes = line.split()

            # Extract the information
            pdb_code = columes[0]
            ligand_name = columes[7][1:-1] # skip ()
            pdbbind_index[pdb_code] = ligand_name
    
    folder_path = "/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed"
    # pdb_codes = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]
    pdb_codes = []
    # df_last = pd.read_csv('/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/Qbio_dataset.csv')
    # df_failed = df_last[df_last['Gen'] == 0]
    # df_success = df_last[df_last['Gen'] == 1]
    
    # to_drop = []

    # for idx, row in tqdm(df_success.iterrows()):
    #     try:
    #         if not os.path.isfile(row['Ligand_file']):
    #             to_drop.append(idx)
    #             continue
    #         if calculate_distance(row['Ligand_file'], row['PDB_file']):
    #             to_drop.append(idx)
    #             os.remove(row['Ligand_file'])
    #     except:
    #         to_drop.append(idx)
    #         os.remove(row['Ligand_file'])
    
    # df_success.drop(to_drop)
    # df_all = df_success.append(df_failed)
    
    # df_all.to_csv('Qbio_dataset_aft.csv', index=False)
    # dic_dif = {}
    # for idx, row in df_all.iterrows():
    #     if row['PDB'] not in dic_dif:
    #         dic_dif[row['PDB']] = 1
    #     else:
    #         dic_dif[row['PDB']] += 1
            
    # total_count = len(dic_dif)
    # greater_than_one_count = sum(value > 1 for value in dic_dif.values())
    # proportion = greater_than_one_count / total_count

    # print(f"The proportion of dic_dif values greater than 1 is: {proportion}")


    # for idx, row in df_failed.iterrows():
    #     if os.path.isfile(row['Ligand_file']):
    #         os.remove(row['Ligand_file'])
    #     if os.path.isfile(row['PDB_file']):
    #         os.remove(row['PDB_file'])
    #     pdb_codes.append(row['PDB'])

    # Some complexes were removed from PDB
    removed = {'4jdf':'7oyz', '4pox':'5oyz', '1qon': '6xyu', '5ab1': '6yhw'}
    removed_back = {'7oyz':'4jdf', '5oyz':'4pox', '6xyu': '1qon', '6yhw': '5ab1'}

    if '1a50' in pdb_codes:
        pdb_codes.remove('1a50') # for some reason, this directory is empty


    print("Submitting tasks to the thread pool...")

    df = pd.read_csv('/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/Q-BioLiP_all.csv', low_memory=False)
    df = df[df['PDB ID'].isin(pdb_codes)]
    multi_bs = []
    # pdb_codes = ['1uys']

    for pdbcode in tqdm(pdb_codes):
        if pdbcode in removed_back:
            pdbcode = removed_back[pdbcode]
        ligand_name = pdbbind_index[pdbcode]

        # in Qbiolib, "-mer" is renamed to "III" ---> Peptide or "k-mer"
        # refer to ligand.json
        if ligand_name[-4:] == "-mer":
            ligand_name = "k-mer"
        handle(pdbcode, ligand_name)
        
        # executor.submit(handle, pdbcode, ligand_name)

    executor.shutdown(wait=True)
    # print(df_result)
    # df_result  = df_result.append(df_success)
    df_result.to_csv('Qbio_dataset.csv', index=False)
    print("FINISHED")
    # print(multi_bs)
    with open('multi_bs.txt', 'w') as f:
        for item in multi_bs:
            f.write("%s\n" % item)