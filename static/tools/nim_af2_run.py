#!/usr/bin/env python3
import os, sys
import requests
import time
from pathlib import Path
import pandas as pd
import numpy as np


def main(seq, name, result_folder, key):
    # Variables
    url = os.getenv("URL", "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2")
    status_url = os.getenv("STATUS_URL", "https://health.api.nvidia.com/v1/status")

    sequence = (seq)
    output_file = Path(f"{result_folder}/{name}.json")

    # Initial request
    headers = {
        "content-type": "application/json",
        "Authorization": f"Bearer {key}",
        "NVCF-POLL-SECONDS": "300",
    }
    data = {
        "sequence": sequence,
        "algorithm": "mmseqs2",
        "e_value": 0.0001,
        "iterations": 1,
        "databases": ["uniref90", "small_bfd", "mgnify"],
        "relax_prediction": False,
        "skip_template_search" : True
    }

    print("Making request...")
    response = requests.post(url, headers=headers, json=data)

    # Check the status code
    if response.status_code == 200:
        output_file.write_text(response.text)
        print(f"Response output to file: {output_file}")
    elif response.status_code == 202:
        print("Request accepted...")
        # Extract reqId header
        req_id = response.headers.get("nvcf-reqid")

        # Poll the /status endpoint
        while True:
            print("Polling for response...")
            status_response = requests.get(f"{status_url}/{req_id}", headers=headers)

            if status_response.status_code != 202:
                output_file.write_text(status_response.text)
                print(f"Response output to file: {output_file}")
                break
    else:
        print(f"Unexpected HTTP status: {response.status_code}")
        print(f"Response: {response.text}")

def getKey():
    origin_path = Path(__file__).resolve().parent.parent.parent
    with open(f"{origin_path}/key.txt", "r") as f:
        key = f.read()
    return key


if __name__ == "__main__":
    result_folder = "/opt/platform/smart_bench/static/output/nim/af3"
    csv_file = "/opt/platform/smart_bench/static/tools/af2_sample.csv"
    
    df = pd.read_csv(csv_file)
    
    for i in range(len(df)):
        name = df.loc[i, "name"]
        seq = df.loc[i, "seqeunce"]
        key = getKey()
        main(seq, name, result_folder, key)
    
    

