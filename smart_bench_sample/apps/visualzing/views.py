import requests
from django.shortcuts import render
from django.core.files.storage import default_storage
from django.conf import settings
from django.http import JsonResponse

from apps.visualzing.services import visualzingFolderGenerator, runPlip, run2DImg, runGninaSocre

import os, json
import pandas as pd
import subprocess
from datetime import datetime
import numpy as np

from apps import getKey, outputPath

# Create your views here.

def visualzingInputPage(request):
    return render(request, 'visualzing/input.html')

def visualzingProcess(request):
    if request.method == 'POST':
        current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
        main_results_dir = f"{outputPath()}/visualzing/{current_time}"
        os.makedirs(main_results_dir, exist_ok=True)

        pdb_file = request.FILES.get('pdb_file')
        pdb_id = request.POST.get('pdb_id')

        # 파일 이름 결정
        if pdb_file:
            filename = pdb_file.name  # 원래 파일 이름
        elif pdb_id:
            filename = f"{pdb_id.upper()}.pdb"  # 입력한 PDB ID 기반 이름
        else:
            return JsonResponse({'error': "No valid PDB input provided."})

        main_file_path = os.path.join(main_results_dir, filename)
        print(main_file_path)
        main_file_name, _ = os.path.splitext(os.path.basename(main_file_path))

        # 파일 저장
        if pdb_file:
            with open(main_file_path, 'wb') as f:
                for chunk in pdb_file.chunks():
                    f.write(chunk)

        elif pdb_id:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = requests.get(url)
            if response.status_code == 200:
                with open(main_file_path, 'wb') as f:
                    f.write(response.content)
            else:
                return JsonResponse({'error': f"PDB ID {pdb_id} not found in RCSB."})
        
        plip_result_dir, img2d_result_dir, gnina_result_dir = visualzingFolderGenerator(main_results_dir)
        
        runPlip(main_file_path, plip_result_dir)
        run2DImg(main_file_path, main_file_name, img2d_result_dir)
        gnina_score_csv_path = runGninaSocre(main_file_path, main_file_name, gnina_result_dir)

        return JsonResponse({'mainResultsDir': main_results_dir})
        
def visualzingResultPage(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    pdb_file = [
        f for f in os.listdir(main_results_dir)
        if f.lower().endswith('.pdb')
    ][0]
    
    with open(f"{main_results_dir}/{pdb_file}", 'r') as f1:
        origin_pdb = f1.read()
    return render(request, 'visualzing/result.html', {"origin_pdb" : origin_pdb})

def visualzingPdb(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)

    pdb_file = [
            f for f in os.listdir(main_results_dir)
            if f.lower().endswith('.pdb')
        ][0]
    return JsonResponse({"origin_pdb" : pdb_file})

def visualzingResult2DImage(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)

    img2d_result_dir = f"{main_results_dir}/img2d"
    png_file = [
        f for f in os.listdir(img2d_result_dir)
        if f.lower().endswith('.png')
    ][0]
    return JsonResponse({'img2D': png_file})

def visualzingPlipResultImages(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)

    print(main_results_dir)

    plip_result_dir = f"{main_results_dir}/plip"
    png_files = [
        f for f in os.listdir(plip_result_dir)
        if f.lower().endswith('.png')
    ][0]
    return JsonResponse({'plip_image': png_files})

# def visualzingPlipResultPdbs(request):
#     main_results_dir = request.GET.get('mainResultsDir')
#     if not main_results_dir or not os.path.isdir(main_results_dir):
#         return JsonResponse({'error': 'Invalid results_dir'}, status=400)

#     plip_result_dir = f"{main_results_dir}/plip"
#     pdb_files = [
#             f for f in os.listdir(plip_result_dir)
#             if f.lower().endswith('.pdb')
#         ][0]
#     return JsonResponse({"plip_pdb" : pdb_files})

def visualzingPlipScore(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    gnina_result = os.path.join(main_results_dir, "plip", "plip.csv")

    try:
        df = pd.read_csv(gnina_result)
        data = df.replace({np.nan: "-"}).to_dict(orient='records')  # [{col1: val1, col2: val2}, ...]
        columns = list(df.columns)
        return JsonResponse({'columns': columns, 'data': data})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)
    
def visualzingPlipLigand(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    plip_ligand_info = os.path.join(main_results_dir, "plip", "ligand_info.json")

    try:
        with open(plip_ligand_info, 'r') as f:
            ligand_info = json.load(f)
        return JsonResponse({'ligand_info': ligand_info})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)

def visualzingGninaScore(request):
    main_results_dir = request.GET.get('mainResultsDir')
    if not main_results_dir or not os.path.isdir(main_results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)
    
    gnina_result = os.path.join(main_results_dir, "gnina", "gnina_score.csv")

    try:
        df = pd.read_csv(gnina_result)
        data = df.to_dict(orient='records')  # [{col1: val1, col2: val2}, ...]
        columns = list(df.columns)
        return JsonResponse({'columns': columns, 'data': data})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)
    
    