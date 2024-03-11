import requests
import os
import sys

import threading
lock = threading.Lock()

import concurrent.futures

arr = []
lock2 = threading.Lock()
k = 0

def version(pdb_code):
    global arr
    url = f"https://www.rcsb.org/versions/{pdb_code}"
    response = requests.get(url)
    if response.status_code != 200:
        return
    response_text = response.text

    # Split the response by "><" to find the date
    split_response = response_text.split('><')
    flag = False
    # Find the first date in the format %d%d%d%d-%d%d-%d%d from the end
    for i in range(len(split_response) - 1, -1, -1):
        date = split_response[i]
        date = date[-14:-4]
        # print(len(date))
        if len(date) == 10 and date[4] == '-' and date[7] == '-':
            flag = True
            break
    if flag and date > '2020-01-01': # Find matched one and date is after 2020-01-01
        with lock:
            arr.append(pdb_code)
    with lock2:
        global k
        k += 1
        print(f"{k}/{len(pdb_codes)}")
        # print(date)  # Replace with your desired logic

    

if __name__ == '__main__':
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=(os.cpu_count()))
    folder_path = "/home/t-kaiyuangao/Workspace/DynamicBind/PDBBindData/PDBBind_processed"
    pdb_codes = [name for name in os.listdir(folder_path) if os.path.isdir(os.path.join(folder_path, name))]

    for i in pdb_codes:
        executor.submit(version, i)
    executor.shutdown(wait=True)

    print(len(arr))
    with open('versionchanged.txt', 'w') as f:
        for i in arr:
            f.write(i+'\n')