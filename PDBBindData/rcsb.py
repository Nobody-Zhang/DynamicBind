import re
import requests
import pandas as pd
import os
import threading

from gemmi import cif
import shutil
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)

import concurrent.futures

from tqdm import tqdm

cnt = 0
success_cnt = 0
tot = 0

lock = threading.Lock()
lock_cnt = threading.Lock()
lock_ass = threading.Lock()
ass_list = []
df = pd.DataFrame()

def find_assembly_index(pdb_code):
    # Read the cif file
    cif_path = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/.tmp/{pdb_code}.cif'
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
        with open(f"/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_code}/{pdb_code}_protein_processed.pdb", "r") as f:
            # Read the first line which contains chains selected
            lines = f.readlines()
            chains = [line.split()[4] for line in lines if line.startswith('ATOM')]
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
    if not flag and len(list_assembly_all) > 1:
        print(f"Warning: {pdb_code} has no matched assembly chains, use the first one")
        with lock_ass:
            global ass_list
            ass_list.append(pdb_code)
    return checked



def download_mol2(pdb_code, ligand_name):
    global df
    global lock
    global mul_cnt 
    global success_cnt
    global cnt
    global failed
    cif_path = f'/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/.tmp/{pdb_code}.cif'
    url_ = 0
    doc = cif.read_file(cif_path)
    block = doc.sole_block()
    # Find assembly Chains
    list_assembly_all = list(block.find_values('_pdbx_struct_assembly_gen.asym_id_list'))
    list_assembly_all = [_.split(',') for _ in list_assembly_all]

    assembly_index = 0
    try:
        assembly_index = find_assembly_index(pdb_code)
    except:
        try:
            assembly_index = find_assembly_index(pdb_code)
        except:
            print(f"Error: {pdb_code} failed to get the assembly chains")
            failed.append(pdb_code)
            return -1
    
    list_assembly = list_assembly_all[assembly_index]
    info_url = f"https://www.rcsb.org/structure/{pdb_code}"
    res_url = requests.get(info_url)
    pattern = f"https://models.rcsb.org/v1/{pdb_code}/ligand\?auth_seq_id=.*_{ligand_name}.mol2"

    split_html = res_url.text.split("><")
    sta = -1

    for _ in split_html:
        matches = re.findall(pattern, _)
        if len(matches) > 0:
            url_ = matches[0].replace("amp;", "")
            parts = url_.split('_')
            chain_ = parts[-2]
            if chain_ in list_assembly:
                dir = f"PDBBind_rcsb/{pdb_code}"
                lig_file = os.path.join(dir, f"{ligand_name}_{chain_}.mol2")
                resp_mol2 = requests.get(url_)
                if resp_mol2.status_code == 200:
                    os.makedirs(dir, exist_ok=True)
                    with open(lig_file, 'w') as f:
                        sta = 1
                        f.write(resp_mol2.text)
                        # copy the pdbbind protein file to the .tmp folder
                        shutil.copy(f"/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_code}/{pdb_code}_protein_processed.pdb", dir)
                else:
                    sta = 0
                with lock:
                    df = df.append({'pdb_code': pdb_code, 'ligand_name': ligand_name, 
                                    'status': sta, 'ligand_path': lig_file,
                                    'protein_path': f'{dir}/{pdb_code}.pdb'}, ignore_index=True)
    
    # cannot find the download url according to ligand name
    if sta == -1: 
        with lock:
            df = df.append({'pdb_code': pdb_code, 'ligand_name': ligand_name, 
                            'status': sta}, ignore_index=True)
    with lock_cnt:
        cnt += 1
        if sta != -1:
            success_cnt += 1
            num_files = len([name for name in os.listdir(f"/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_rcsb/{pdb_code}") if os.path.isfile(os.path.join(f"PDBBind_rcsb/{pdb_code}", name))])
            if num_files > 2:
                mul_cnt += 1
        if success_cnt != 0:
            print(f"MULTIPLE MOL PERCENTAGE: {mul_cnt/success_cnt*100}%")

        print(f"Finished: {cnt}, Success: {success_cnt}, Failed: {cnt - success_cnt}, Total: {tot}, Percentage: {cnt/tot*100}%")
        
    return 0

if __name__ == '__main__':
    # Thread pool
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=(os.cpu_count()))
    # dictionary to store the mapping information from pdb_code to ligand_name, extracted from INDEX_general_PL_data.2020
    pdbbind_index = {}
    # Parse PDBBind index
    with open('INDEX_general_PL_data.2020', 'r') as f:
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
    pdb_codes = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]
    
    # rscb has removed (renamed) some complexes, we need to remove them from the list
    # we need to rename the pdbcode to the newer one
    # !! TASKS TO BE FINISHED !!
    removed = []

    removed.append('4jdf')
    removed.append('4pox')
    removed.append('1qon')
    removed.append('5ab1')
    

    pdb_codes.remove('1a50') # for some reason, this directory is empty

    failed = []
    mul_cnt = 0
    for pdb_code in tqdm(pdb_codes):
        # There are some complexes removed from rscb. Mannual check needed.
        if pdb_code in removed: 
            continue
        
        tot += 1
        ligand_name = pdbbind_index[pdb_code]
        executor.submit(download_mol2, pdb_code, ligand_name)
        
    print("Now, all tasks submitted")
    print("\n\n=========================\n\n")
    executor.shutdown(wait=True)

    with lock:
        df.to_csv('pdbbind_rcsb.csv', index=False)
        print("csv saved in pdbbind_rcsb.csv")
    
    with open('failed_assembly_chains.txt', 'w') as f:
        for item in failed:
            f.write("%s\n" % item)

    with open('failed_assembly_choose_1.txt', 'w') as f:
        for i in ass_list:
            f.write("%s\n" % i)
    print("All finished")
    