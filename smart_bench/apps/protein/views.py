from django.shortcuts import render
from django.http import JsonResponse
from django.utils.http import urlsafe_base64_decode
import requests
from pathlib import Path
import json, os
import pandas as pd
import numpy as np
import tempfile
from datetime import datetime
import subprocess
import time
import shutil
from urllib.parse import unquote
from concurrent.futures import ProcessPoolExecutor, as_completed

from apps import getKey, outputPath
from apps.protein.services import setSample, complexSample, runDownloadPDB, proteinMakeResultFolder, runOpenfold, runAlphafold2, runEsmfold, runBoltzChai, proteinAlignment
# Create your views here.

def proteinInputPage(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    samples = setSample()
    # print(samples)
    return render(request, 'protein/input.html', {'databases': databases, 'samples_json': json.dumps(samples)})

def proteinInputPageSample(request):
    databases = ['Uniref30_2302', 'PDB70_220313', 'colabfold_envdb_202108']
    samples = setSample()
    # print(samples)
    return render(request, 'protein/input_sample.html', {'databases': databases, 'samples_json': json.dumps(samples)})

def proteinProcess(request):
    if request.method == "POST":    
        try:
            current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
            results_dir = f"{outputPath()}/protein/{current_time}"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')
            protein = data.get('aminoSeq', '')
            msa_db = data.get('msaDatabases', '')
            print("Input protein sequence: ", protein)
            
            os.makedirs(results_dir, exist_ok=True)
            openfold_result_dir, alphafold_result_dir, esmfold_result_dir, boltz_result_dir, chai_result_dir = proteinMakeResultFolder(results_dir)
            # print(opnefold_result)
            
            # runOpenfold(openfold_result_dir, msa_db, protein)
            # runAlphafold2(alphafold_result_dir, protein)
            # runEsmfold(esmfold_result_dir, protein)
            # runBoltzChai(results_dir, boltz_result_dir, chai_result_dir, protein)
            
            print("Origin protein PDB download.....")
            runDownloadPDB(results_dir, sample_name)
            
            # 병렬 실행
            with ProcessPoolExecutor(max_workers=4) as executor:
                futures = {
                    executor.submit(runOpenfold, openfold_result_dir, msa_db, protein): "Openfold",
                    executor.submit(runAlphafold2, alphafold_result_dir, protein, sample_name): "Alphafold2",
                    executor.submit(runEsmfold, esmfold_result_dir, protein): "Esmfold",
                    executor.submit(runBoltzChai, results_dir, boltz_result_dir, chai_result_dir, protein): "Boltz&Chai"
                }

                for future in as_completed(futures):
                    task_name = futures[future]
                    try:
                        result = future.result()
                        print(f"{task_name} 완료: {result}")
                    except Exception as e:
                        print(f"{task_name} 실패: {e}")
                        raise  # 예외를 다시 발생시켜서 전체 중단
            
            Structure_alignment = proteinAlignment(results_dir)
            print(Structure_alignment)

            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def proteinProcessSample(request):
    if request.method == "POST":    
        try:
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')
            results_dir = f"{outputPath()}/protein/_{sample_name}"

            return JsonResponse({"status": "done", "results_dir" : results_dir})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def proteinOutputPage(request):
    results_dir = request.GET.get('results_dir')
        
    # with open(f"{results_dir}/rot-origin.pdb", 'r') as f1:
    #     origin = f1.read()
        
    pdb_map = {
        "origin_pdb": "rot-origin.pdb",
        "openfold_pdb": "rot-openfold.pdb",
        "af2_pdb": "rot-alphafold.pdb",
        "esmfold_pdb": "rot-esmfold.pdb",
        "boltz_pdb": "rot-boltz.pdb",
        "chai_pdb": "rot-chai.pdb",
    }

    context = {
        key: open(f"{results_dir}/{filename}", 'r').read()
        for key, filename in pdb_map.items()
    }

    return render(request, 'protein/result.html', context)

def boltzInputPage(request):
    samples = setSample()
    # print(samples)
    return render(request, 'boltz/input.html', {'samples_json': json.dumps(samples)})

def boltzProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltzOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print(sample_name)
    sample_pdb_path = os.path.join(results_dir, f"{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'boltz/result.html', {
        "result_pdbs": pdb_content
    })
    
def chaiInputPage(request):
    samples = setSample()
    # print(samples)
    return render(request, 'chai/input.html', {'samples_json': json.dumps(samples)})

def chaiProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/chai"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def chaiOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    print(sample_name)
    sample_pdb_path = os.path.join(results_dir, f"{sample_name}.pdb")

    pdb_content = ""
    if os.path.exists(sample_pdb_path):
        with open(sample_pdb_path, 'r') as f:
            pdb_content = f.read()

    return render(request, 'chai/result.html', {
        "result_pdbs": pdb_content
    })
    
# boltz complex 
def boltzComplexInputPage(request):
    samples = complexSample()
    # print(samples)
    return render(request, 'boltz/complex_input.html', {'samples_json': json.dumps(samples)})

def boltzComplexProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/boltz_complex"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            print(sample_name)
            time.sleep(1.5)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def boltzComplexOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')
    
    # print(sample_name)
    sample_json_path = os.path.join(results_dir, f"{sample_name}.json")

    pdb_content = ""
    if os.path.exists(sample_json_path):
        with open(sample_json_path, 'r') as f:
            pdb_content = json.load(f)

    return render(request, 'boltz/complex_result.html', {
        "result_pdbs": pdb_content
    })
    

def chaiComplexInputPage(request):
    samples = complexSample()
    # print(samples)
    return render(request, 'chai/complex_input.html', {'samples_json': json.dumps(samples)})

def chaiComplexProcess(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/chai_complex"
            
            data = json.loads(request.body)
            sample_name = data.get('sampleName', '')    
            time.sleep(1.5)  # 5초 동안 멈춤
            return JsonResponse({"status": "done", "results_dir" : results_dir, "sample_name" : sample_name})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def chaiComplexOutputPage(request):
    sample_name = request.GET.get('sample_name')
    results_dir = request.GET.get('results_dir')

    # print(sample_name)
    sample_json_path = os.path.join(results_dir, f"{sample_name}.json")

    pdb_content = ""
    if os.path.exists(sample_json_path):
        with open(sample_json_path, 'r') as f:
            pdb_content = json.load(f)

    return render(request, 'chai/complex_result.html', {
        "result_pdbs": pdb_content
    })