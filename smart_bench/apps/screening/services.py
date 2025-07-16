import requests
import json
from pathlib import Path
import json, os
import subprocess

# def nvidiaOpenfold(protein):
#     OPENFOLD2_HOST = 'http://lnplicense.iptime.org:8701'
#     of2_response = requests.post(
#         f'{OPENFOLD2_HOST}/biology/openfold/openfold2/predict-structure-from-msa-and-template',
#         json={
#             'sequence': protein,
#             'use_templates': False,
#             'relaxed_prediction': False,
#             # 'alignments': msa_response['alignments'],
#     }).json()

#     folded_protein = of2_response["structures_in_ranked_order"].pop(0)["structure"]
#     return folded_protein

def setChemSample():
    samples = [
        {"name": "SODC", "smiles": "CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"},
        {"name": "COX2", "smiles": "CCOC(=O)[C@H](CCC1=CC=CC=C1)N[C@@H](C)C(=O)N2CC3(C[C@H]2C(=O)O)SCCS3"},
        {"name": "TNFA", "smiles": "CC(C)NC1=NC2=CC(=C(C=C2N1[C@@H]3[C@H]([C@H]([C@@H](O3)CO)O)O)Cl)Cl"},
    ]
    return samples

def nvidiaOpenfold(protein, msa, results_dir):
    key = "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"
    url = os.getenv("URL", "https://health.api.nvidia.com/v1/biology/openfold/openfold2/predict-structure-from-msa-and-template")
    output_file = Path(f"{results_dir}/openfold.json")
    selected_models = [1]
    sequence = (protein)

    data = {
        "sequence": sequence,
        "alignments": msa['alignments'],
        "selected_models": selected_models,
        "relax_prediction": False,
    }
    print(data)

    # ---------------------------------------------------------
    # Submit
    # ---------------------------------------------------------
    headers = {
        "content-type": "application/json",
        "Authorization": f"Bearer {key}",
        "NVCF-POLL-SECONDS": "300",
    }
    print("Making request...")
    response = requests.post(url, headers=headers, json=data)
    print(response)

    # ---------------------------------------------------------
    # View response
    # ---------------------------------------------------------
    if response.status_code == 200:
        output_file.write_text(response.text)
        print(f"Response output to file: {output_file}")
        with open(f"{results_dir}/openfold.json", 'r') as f:
            openfold = json.load(f)
            folded_protein = openfold['structures_in_ranked_order'][0]['structure']
            print(openfold)
    else:
        print(f"Unexpected HTTP status: {response.status_code}")
        print(f"Response: {response.text}")
        
    return folded_protein

def nvidiaMsa(MSA_HOST):
    pass

def nvidiaMolmin(smiles, num_mols, property_name, minimize, min_similarity, particles, iterations):
    invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/molmim/generate"

    headers = {
        "Authorization": "Bearer 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",
        "Accept": "application/json",
    }

    payload = {
    "algorithm": "CMA-ES",
    "num_molecules": num_mols,
    "property_name": property_name,
    "minimize": minimize,
    "min_similarity": min_similarity,
    "particles": particles,
    "iterations": iterations,
    "smi": smiles
    }

    # re-use connections
    session = requests.Session()
    print(session)
    _molmin_response = session.post(invoke_url, headers=headers, json=payload)
    print(_molmin_response)
    _molmin_response.raise_for_status()
    print(_molmin_response)
    molmin_response = json.loads(_molmin_response.json()['molecules'])
    print(_molmin_response.json())
    return molmin_response

# def nvidiaMolmin(smiles, num_mols, property_name, minimize, min_similarity, particles, iterations):
#     MOLMIN_HOST = 'http://lnplicense.iptime.org:8705'
    
#     # print(type(MOLMIN_HOST), type(smiles), type(num_mols), type(property_name), type(minimize), type(min_similarity), type(particles), type(iterations))
#     headers = {
#         "Accept": "application/json",
#         "Content-Type": "application/json",
#     }
    
#     data = {
#         "algorithm": "CMA-ES",
#         "num_molecules": num_mols,
#         "property_name": property_name,
#         "minimize": minimize,
#         "min_similarity": min_similarity,
#         "particles": particles,
#         "iterations": iterations,
#         "smi": smiles
#     }
    
#     response = requests.post(f"{MOLMIN_HOST}/generate", headers=headers, data=json.dumps(data))
#     molmin_response = response.json()
#     # smiles_list = '\n'.join([ v['smiles']for v in molmin_response['generated']])
#     return molmin_response

def nvidiaDiffdock(folded_protein, generated_ligands, poses):
    DIFFDOCK_HOST = 'http://lnplicense.iptime.org:8703'
    diffdock_response = requests.post(
    f'{DIFFDOCK_HOST}/molecular-docking/diffdock/generate',
    json={
        'protein': folded_protein,
        'ligand': generated_ligands,
        'ligand_file_type': 'txt',
        'num_poses': poses,
        'time_divisions': 20,
        'num_steps': 18,
    }).json()
    
    return diffdock_response
