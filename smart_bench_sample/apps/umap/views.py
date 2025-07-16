from django.shortcuts import render
from django.http import JsonResponse

import time
from datetime import datetime
import json, os
import pandas as pd
import numpy as np
import tempfile
import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
from django.http import HttpResponse

from apps import getKey, outputPath

from apps.umap.services import molMin, genMol, molsparkProcess
# Create your views here.

def umapInputPage(request):
    sections = ["MolMin", "GenMol", "MolSpark"]
    return render(request, 'umap/input.html', {'sections': sections})

def umapInputPageSample(request):
    sections = ["MolMin", "GenMol", "MolSpark"]
    return render(request, 'umap/input_sample.html', {'sections': sections})

def generateMolecules(request):
    if request.method == "POST":    
        try:
            current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
            results_dir = f"{outputPath()}/umap/{current_time}"
            
            data = json.loads(request.body)  
            smiles = data.get('molSeq', '')
            gen_smiles = ""

            if smiles == "CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3":
                gen_smiles = "C124CN3C1.S3(=O)(=O)CC.C4C#N.[*{20-20}]"
            elif smiles == "CCOC(=O)[C@H](CCC1=CC=CC=C1)N[C@@H](C)C(=O)N2CC3(C[C@H]2C(=O)O)SCCS3":
                gen_smiles = "N13CC2(CC14)SCCS2.C4(=O)O.[*{20-25}]"
            else:
                gen_smiles = "C12OC3C(O)C1O.C3O.[*{25-25}]"
            
            # 1. 입력한 smiles를 데이터프레임으로 변환
            input_df = pd.DataFrame([{"smiles": smiles, "tool" : "input"}])
            
            # 함수 호출
            molmin_df = molMin(smiles)
            genmol_df = genMol(gen_smiles)
            molspark_df = molsparkProcess(smiles)
            
            os.makedirs(results_dir, exist_ok=True)
            result_df = pd.concat([input_df, genmol_df, molmin_df, molspark_df], axis=0)
            result_df.to_csv(f"{results_dir}/gen3tools.csv", index=False)
            _, result_df_csv = tempfile.mkstemp(suffix=".csv")
            result_df.to_csv(result_df_csv, index=False)
            
            subprocess.run([f"/opt/anaconda3/envs/chemplot/bin/python /opt/platform/smart_bench/static/tools/computation_property.py {result_df_csv} {results_dir}"], shell=True)
            # subprocess.run([f"/opt/anaconda3/envs/chemplot/bin/python /opt/platform/smart_bench/static/tools/computation_property.py {results_dir}/gen3tools.csv {results_dir}"], shell=True)
            
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def generateMoleculesSample(request):
    if request.method == "POST":    
        try:
            current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
            sample_dir = f"{outputPath()}/umap_sample"
            
            data = json.loads(request.body)  
            sample = data.get('molSeq', '')
            results_dir = f"{sample_dir}/{sample}"
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def generate3dSdf(request):
    smiles = request.GET.get('smiles')
    if not smiles:
        return HttpResponse("Missing SMILES", status=400)

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return HttpResponse("Invalid SMILES", status=400)

    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

    sdf = Chem.MolToMolBlock(mol)
    return HttpResponse(sdf, content_type='chemical/x-mdl-sdfile')
       
def umapResultPage(request):
    umap_result = request.GET.get('umap_result')
    if not umap_result or not os.path.isdir(umap_result):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    result_cal_df = pd.read_csv(f"{umap_result}/Rdkit_calculate_result.csv").drop(columns=[f"fp_{i+1}" for i in range(2048)])
    return render(request, 'umap/result.html',  {"smiles" : result_cal_df['smiles'].tolist(), "cal_data": result_cal_df.replace({np.nan: None}).to_dict(orient="records")})

def umapResultToolNum(request):
    umap_result = request.GET.get('umap_result')
    if not umap_result or not os.path.isdir(umap_result):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    result_df = pd.read_csv(f"{umap_result}/gen3tools.csv")
    # os.remove(f"{umap_result}/gen3tools.csv")
    genmol_num, molmin_num, molspark_num = len(result_df[result_df.loc[:, "tool"] == "GenMol"]), len(result_df[result_df.loc[:, "tool"] == "MolMin"]), len(result_df[result_df.loc[:, "tool"] == "MolSpark"])
    return JsonResponse({"genmol_num": genmol_num, "molmin_num" : molmin_num, "molspark_num" : molspark_num}) 

def umapResultImage(request):
    umap_result = request.GET.get('umap_result')
    if not umap_result or not os.path.isdir(umap_result):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    umap_img = [
        f for f in os.listdir(umap_result)
        if f.lower().endswith('.png')
    ]
    return JsonResponse({'umap_img': umap_img})

def umapResultCsv(request):
    umap_result = request.GET.get('umap_result')
    if not umap_result or not os.path.isdir(umap_result):
        return JsonResponse({'status': 'error', 'message': 'Invalid results_dir'}, status=400)

    try:
        filename = "Rdkit_calculate_result.csv"
        full_path = os.path.join(umap_result, filename)

        if not os.path.isfile(full_path):
            return JsonResponse({'status': 'error', 'message': 'File not found'}, status=404)

        # 프론트에서 static 기준으로 접근할 수 있게 상대 경로 전달
        relative_path = full_path.split("/static/")[-1]  # static부터 시작
        download_url = f"/static/{relative_path}"

        return JsonResponse({
            "status": "success",
            "download_url": download_url
        })
    except Exception as e:
        return JsonResponse({
            "status": "error",
            "message": str(e)
        })

# def umapResultCsv(request):
#     umap_result = request.GET.get('umap_result')
#     if not umap_result or not os.path.isdir(umap_result):
#         return JsonResponse({'error': 'Invalid results_dir'}, status=400)

#     try:
#         df = pd.read_csv(f"{umap_result}/cal_result.csv")

#         # 필요한 경우 데이터 전처리
#         df.fillna('-', inplace=True)  # null 처리

#         data = df.to_dict(orient='records')  # [{...}, {...}] 형태
#         return JsonResponse({
#             "status": "success",
#             "data": data
#         })
#     except Exception as e:
#         return JsonResponse({
#             "status": "error",
#             "message": str(e)
#         })