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

cnt = 0
success_cnt = 0
tot = 0

lock = threading.Lock()
lock_cnt = threading.Lock()
df = pd.DataFrame()

def test_cp(pdb_code, test_mode=False):
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
            chains = lines[0].split("'")[1].replace("chain ", "").split(" or ")
            chains_tmp = []
            for idx, chain in enumerate(chains):
                if len(chain) > 1: # This means that there are string like ... or

                    las = chain[-1]
                    fir = chains[idx - 1][0]
                    
                    
                    for t in range(ord(fir) + 1, ord(las) + 1):
                        chains_tmp.append(chr(t))
                else:
                    chains_tmp.append(chain)
            chains = chains_tmp
            if test_mode:
                print("chains: ", chains)
                print("dic_name: ", dic_name)
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
    # checked = 0
    # flag = False
    # if len(list_assembly_all) > 1:
    #     # Means that there are multiple assemblies,
    #     # Need information from the pdbbind website
    #     sequence_url = f"https://www.rcsb.org/annotations/{pdb_code}"
    #     resp_seq = requests.get(sequence_url)
    #     # The annotations website guarantees that there is only protein information
    #     pattern_auth = re.compile(r'\[auth')
    #     matches_auth = pattern_auth.finditer(resp_seq.text)

    #     positions = [match.start() for match in matches_auth]
    #     dic_name = {}
    #     if positions != []:
    #         # This means that there are changes of the chain name
    #         for _ in positions:
    #             pdb_name = resp_seq.text[_+7]
    #             cif_name = resp_seq.text[_-2]
    #             dic_name[pdb_name] = cif_name
    #     with open(f"PDBBindData/PDBBind_processed/{pdb_code}/{pdb_code}_protein_processed.pdb", "r") as f:
    #         # Read the first line which contains chains selected
    #         lines = f.readlines()
    #         chains = lines[0].split("'")[1].replace("chain ", "").split(" or ")
    #         chains_tmp = []
    #         for idx, chain in enumerate(chains):
    #             if len(chain) > 1: # This means that there are string like ... or

    #                 las = chain[-1]
    #                 fir = chains[idx - 1][0]
                    
                    
    #                 for t in range(ord(fir) + 1, ord(las) + 1):
    #                     chains_tmp.append(chr(t))
    #             else:
    #                 chains_tmp.append(chain)
    #         chains = chains_tmp
    #         # Now we know the chains contained in the pdb file
    #         for idx, chain in enumerate(chains):
    #             if chain in dic_name:
    #                 chains[idx] = dic_name[chain]
    #         for idx, chain_list in enumerate(list_assembly_all):
    #             flag = True
    #             for chain in chains:
    #                 if chain not in chain_list:
    #                     flag = False
    #                     break
    #             if flag:
    #                 # Find the matched chain
    #                 checked = idx
    #                 break
    # if not flag and len(list_assembly_all) > 1:
    #     print(f"Warning: {pdb_code} has no matched assembly chains, use the first one")
    checked = 0
    try:
        checked = test_cp(pdb_code)
    except:
        try:
            checked = test_cp(pdb_code)
        except:
            print(f"Error: {pdb_code} failed to get the assembly chains")
            failed.append(pdb_code)
            return -1
            exit(-1)
    list_assembly = list_assembly_all[checked]
    info_url = f"https://www.rcsb.org/structure/{pdb_code}"
    res_url = requests.get(info_url)
    pattern = f"https://models.rcsb.org/v1/{pdb_code}/ligand\?auth_seq_id=.*_{ligand_name}.mol2"
    split_html = res_url.text.split("><")
    sta = -1
    # sta -1表示没找到相应的pdb文件， 0表示网络问题， 1表示成功下载
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
                        # 拷贝文件
                        shutil.copy(f"/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed/{pdb_code}/{pdb_code}_protein_processed.pdb", dir)
                else:
                    sta = 0
                with lock:
                    df = df.append({'pdb_code': pdb_code, 'ligand_name': ligand_name, 
                                    'status': sta, 'ligand_path': lig_file,
                                    'protein_path': f'{dir}/{pdb_code}.pdb'}, ignore_index=True)
    
    if sta == -1: # 不能通过ligand name找到下载的url
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
    executor = concurrent.futures.ThreadPoolExecutor()
    # dictionary to store the information
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
    
    # rscb has removed some complexes, we need to remove them from the list 
    removed = []

    removed.append('4jdf')
    removed.append('4pox')
    removed.append('1qon')
    removed.append('5ab1')
    

    pdb_codes.remove('1a50') # for some reason, this directory is empty

    failed = []
    tot = len(pdb_codes)
    mul_cnt = 0
    from tqdm import tqdm
    for pdb_code in tqdm(pdb_codes):
        if pdb_code in removed: # There are some complexes removed from rscb. Mannual check needed.
            continue
        
        # tot += 1
        ligand_name = pdbbind_index[pdb_code]
        download_mol2(pdb_code, ligand_name)
    with lock:
        df.to_csv('pdbbind_rcsb.csv', index=False)
        print("csv saved in pdbbind_rcsb.csv")
    print("All finished")
    print(f"Failed: {failed}")

    """
    for pdb_code in pdb_codes:
        if pdb_code in removed: # There are some complexes removed from rscb. Mannual check needed.
            continue
        
        tot += 1
        ligand_name = pdbbind_index[pdb_code]
        executor.submit(download_mol2, pdb_code, ligand_name)
        # download_mol2(pdb_code, ligand_name)
    print("Now, all tasks submitted")
    executor.shutdown(wait=True)
    with lock:
        df.to_csv('pdbbind_rcsb.csv', index=False)
        print("csv saved in pdbbind_rcsb.csv")
    print("All finished")
    """
    
    