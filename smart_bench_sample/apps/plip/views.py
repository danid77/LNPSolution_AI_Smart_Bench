import os
import requests
from django.shortcuts import render
from django.core.files.storage import default_storage
from django.conf import settings
from django.http import JsonResponse

import tempfile
import json
import os
import subprocess
from datetime import datetime
import xml.etree.ElementTree as ET
import pandas as pd

from apps import getKey, outputPath

# Create your views here.

def plipInputPage(request):
    return render(request, 'plip/input.html')

def plipProcess(request):
    if request.method == 'POST':
        current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
        results_dir = f"{outputPath()}/plip/{current_time}/"
        os.makedirs(results_dir, exist_ok=True)

        pdb_file = request.FILES.get('pdb_file')
        pdb_id = request.POST.get('pdb_id')

        # 파일 이름 결정
        if pdb_file:
            filename = pdb_file.name  # 원래 파일 이름
        elif pdb_id:
            filename = f"{pdb_id.upper()}.pdb"  # 입력한 PDB ID 기반 이름
        else:
            return JsonResponse({'error': "No valid PDB input provided."})

        temp_path = os.path.join(results_dir, filename)

        # 파일 저장
        if pdb_file:
            with open(temp_path, 'wb') as f:
                for chunk in pdb_file.chunks():
                    f.write(chunk)

        elif pdb_id:
            url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
            response = requests.get(url)
            if response.status_code == 200:
                with open(temp_path, 'wb') as f:
                    f.write(response.content)
            else:
                return JsonResponse({'error': f"PDB ID {pdb_id} not found in RCSB."})

        # Docker로 PLIP 실행
        try:
            run_plip = subprocess.run(
                ["docker", "run", "--rm",
                 "-v", f"{results_dir}:/results",
                 "-w", "/results",
                 "pharmai/plip:latest", "-f", f"/results/{filename}", "-pxty"],
                capture_output=True, text=True, check=True
            )

        except subprocess.CalledProcessError as e:
            return JsonResponse({
                'error': 'PLIP 분석 중 오류 발생',
                'details': e.stderr
            })

        return JsonResponse({'results_dir': results_dir})


def plipResultPage(request):
    return render(request, 'plip/result.html')

def plipResultImages(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir or not os.path.isdir(results_dir):
        return JsonResponse({'error': 'Invalid results_dir'}, status=400)

    png_files = [
        f for f in os.listdir(results_dir)
        if f.lower().endswith('.png')
    ]
    return JsonResponse({'images': png_files})

def plipResultXml(request):
    results_dir = request.GET.get('results_dir')
    if not results_dir:
        return JsonResponse({'error': 'results_dir 파라미터가 필요합니다.'}, status=400)

    xml_path = os.path.join(results_dir, "report.xml")
    if not os.path.exists(xml_path):
        return JsonResponse({'error': 'XML 파일이 존재하지 않습니다.'}, status=404)

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        sites = []

        for site in root.findall('.//bindingsite'):
            site_data = {
                'id': site.get('id'),
                'ligand': {
                    'longname': site.findtext('./identifiers/longname'),
                    'hetid': site.findtext('./identifiers/hetid'),
                    'chain': site.findtext('./identifiers/chain'),
                    'position': site.findtext('./identifiers/position'),
                    'smiles': site.findtext('./identifiers/smiles'),
                    'inchikey': site.findtext('./identifiers/inchikey'),
                },
                'bs_residues': [],
                'metal_complexes': [],
                'hydrogen_bonds': [],
                'pi_stacks': [],
            }

            for residue in site.findall('.//bs_residues/bs_residue'):
                site_data['bs_residues'].append({
                    'aa': residue.attrib.get('aa'),
                    'id': residue.attrib.get('id'),
                    'contact': residue.attrib.get('contact'),
                    'min_dist': residue.attrib.get('min_dist'),
                    'text': residue.text,
                })

            for metal in site.findall('.//metal_complexes/metal_complex'):
                site_data['metal_complexes'].append({
                    'resnr': metal.findtext('resnr'),
                    'restype': metal.findtext('restype'),
                    'reschain': metal.findtext('reschain'),
                    'metal_type': metal.findtext('metal_type'),
                    'dist': metal.findtext('dist'),
                    'location': metal.findtext('location'),
                    'geometry': metal.findtext('geometry'),
                })

            for hbond in site.findall('.//hydrogen_bonds/hydrogen_bond'):
                site_data['hydrogen_bonds'].append({
                    'resnr': hbond.findtext('resnr'),
                    'restype': hbond.findtext('restype'),
                    'reschain': hbond.findtext('reschain'),
                    'dist_d-a': hbond.findtext('dist_d-a'),
                    'don_angle': hbond.findtext('don_angle'),
                })

            for pi in site.findall('.//pi_stacks/pi_stack'):
                site_data['pi_stacks'].append({
                    'resnr': pi.findtext('resnr'),
                    'restype': pi.findtext('restype'),
                    'reschain': pi.findtext('reschain'),
                    'type': pi.findtext('type'),
                    'centdist': pi.findtext('centdist'),
                    'angle': pi.findtext('angle'),
                })

            sites.append(site_data)

        return JsonResponse({'sites': sites})
    except Exception as e:
        return JsonResponse({'error': str(e)}, status=500)