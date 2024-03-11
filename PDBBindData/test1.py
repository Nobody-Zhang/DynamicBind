from gemmi import cif
import threading
import rcsb
import concurrent.futures
from tqdm import tqdm
import os
import requests

folder_path = "PDBBind_processed"
folder_names = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]

def f(folder_name):
    cif_url = f"https://www.rcsb.org/annotations/{folder_name}"
    # res_url = requests.get(cif_url)
    # matches = re.findall(desir, res_url.text)
    # if len(matches) == 0:
    #     print(i)
    # print(matches)
    response = requests.get(cif_url)
    if response.status_code == 200:
        return
    print(f"Not Found {folder_name}")

executor = concurrent.futures.ThreadPoolExecutor()

for i in tqdm(folder_names):
    executor.submit(f, i)

executor.shutdown(wait=True)

# print(cnt, len(folder_names))
# print(idx)
