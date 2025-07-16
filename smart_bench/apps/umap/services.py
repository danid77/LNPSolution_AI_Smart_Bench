import requests
import json
import numpy as np
import pandas as pd
import os, sys
import tempfile
import subprocess

from apps import getKey

def genMol(smiles):
    invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/genmol/generate"

    headers = {
        "Authorization": f"Bearer {getKey()}",
        "Accept": "application/json",
    }

    payload = {
        "smiles": smiles,
        "num_molecules": 30,
        "temperature": 1,
        "noise": 0,
        "step_size":1,
        "scoring": "QED"
    }

    # re-use connections
    session = requests.Session()
    print(session)
    genmol_response = session.post(invoke_url, headers=headers, json=payload)
    print(genmol_response)
    genmol_response.raise_for_status()
    print(genmol_response.json())
    genmol_df = pd.DataFrame(genmol_response.json()['molecules'])[["smiles"]].replace('', np.nan).dropna(subset="smiles")
    genmol_df['tool'] = "GenMol"
    return genmol_df

def molMin(smiles):
    invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"

    headers = {
        "Authorization": f"Bearer {getKey()}",
        "Accept": "application/json",
    }

    payload = {
    "algorithm": "CMA-ES",
    "num_molecules": 30,
    "property_name": "QED",
    "minimize": False,
    "min_similarity": 0.3,
    "particles": 30,
    "iterations": 10,
    "smi": smiles
    }

    # re-use connections
    session = requests.Session()
    print(session)
    molmin_response = session.post(invoke_url, headers=headers, json=payload)
    print(molmin_response)
    molmin_response.raise_for_status()
    molmin = json.loads(molmin_response.json()['molecules'])
    print(molmin)
    molmin_df = pd.DataFrame(molmin).rename(columns={"sample" : "smiles"})[["smiles"]].replace('', np.nan).dropna(subset="smiles")
    molmin_df["tool"] = "MolMin"
    return molmin_df

def molsparkProcess(smiles):
    # print(smiles, model_type, strategy, temperature, num_samples)
    # 임시 파일 생성
    with tempfile.NamedTemporaryFile(mode="w", delete=False) as smiles_file:
        smiles_file.write(smiles)  # str 그대로 쓰기
        smiles_path = smiles_file.name  # 파일 경로 저장

    print(f"임시 파일 경로: {smiles_path}")

    # 파일이 잘 생성됐는지 확인
    with open(smiles_path, "r") as f:
        print(f.read())

    _, output_file = tempfile.mkstemp(suffix=".csv")
    
    # 기본 JSON 데이터 구조
    reinvent_option = {
        "run_type": "sampling",
        "use_cuda": True,
        "parameters": {
            "model_file": f"/opt/git-tools/REINVENT4/models/mol2mol_high_similarity.prior",
            "smiles_file": smiles_path,
            "sample_strategy": "beamsearch",
            "output_file": output_file,
            "num_smiles": 40,
            "unique_molecules": True,
            "randomize_smiles": True
        }
    }

    # 임시 JSON 파일 생성
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as option_json:
        json.dump(reinvent_option, option_json, indent=4)
        option_json_path = option_json.name

    print(f"임시 JSON 파일 경로: {option_json_path}")
    print(output_file)
    
    reinvent = ["/opt/anaconda3/envs/reinvent4/bin/reinvent", "-f", "json", option_json_path]
    process1 = subprocess.run(reinvent, capture_output=True, text=True)
    print("STDOUT:", process1.stdout)
    print("STDERR:", process1.stderr)
    
    mol2mol_df = pd.read_csv(output_file).head(30)
    
    os.remove(smiles_path)
    os.remove(option_json_path)
    os.remove(output_file)
    
    mol2mol_df = mol2mol_df.rename(columns={"SMILES" : "smiles"})[["smiles"]]
    mol2mol_df['tool'] = "MolSpark"
    # csv_list = os.listdir(cal_output_dir)
    # csv_file = [file for file in csv_list if file.find(".csv") !=-1][0]
    # cal_output_file = f"{cal_output_dir}/{csv_file}"
    
    return mol2mol_df