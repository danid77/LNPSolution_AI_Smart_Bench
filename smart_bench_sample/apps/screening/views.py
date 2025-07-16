from django.shortcuts import render
from django.http import JsonResponse
from django.http import HttpResponse
from django.views.decorators.csrf import csrf_exempt
import requests
import json, os
from datetime import datetime
import subprocess
import time

from apps import getKey, outputPath

from apps.screening.services import setChemSample, nvidiaDiffdock, nvidiaMolmin, nvidiaMsa, nvidiaOpenfold
from apps.protein.services import setSample
# Create your views here.

# MSA_HOST = 'http://localhost:8081'
# OPENFOLD2_HOST = 'http://lnplicense.iptime.org:8701'
# GENMOL_HOST = 'http://lnplicense.iptime.org:8702'
DIFFDOCK_HOST = 'http://lnplicense.iptime.org:8703'
# MOLMIN_HOST = 'http://lnplicense.iptime.org:8705'

def complexInputPage(request):
    sections = ["OpenFold2", "MMseqs2 MSA", "MolMin", "DiffDock Inputs"]
    databases = ['Uniref30_2302', 'colabfold_envdb_202108', 'PDB70_220313']
    protein_samples = setSample()
    # smiles_samples = setChemSample()
    return render(request, 'screening/input.html', {
        'sections': sections,
        'databases': databases,
        'protein_json': json.dumps(protein_samples),
        # 'smiles_json': json.dumps(smiles_samples)
    })

def complexInputPageSample(request):
    sections = ["OpenFold2", "MMseqs2 MSA", "MolMin", "DiffDock Inputs"]
    databases = ['Uniref30_2302', 'colabfold_envdb_202108', 'PDB70_220313']
    protein_samples = setSample()
    # smiles_samples = setChemSample()
    return render(request, 'screening/input_sample.html', {
        'sections': sections,
        'databases': databases,
        'protein_json': json.dumps(protein_samples),
        # 'smiles_json': json.dumps(smiles_samples)
    })

def complexOutputPage(request):
    return render(request, 'screening/result.html')

# def generateMolecules(request):
#     with open('/opt/platform/smart_bench/apps/screening/test_nvidia.json', 'r') as f:
#         data_dict = json.load(f)
#     return JsonResponse(data_dict)


def generateMolecules(request):
    if request.method == "POST":    
        try:
            current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
            results_dir = f"{outputPath()}/complex/{current_time}/"
            os.makedirs(results_dir, exist_ok=True)
            
            data = json.loads(request.body)
            protein = data.get('aminoSeq', '')
            msa_db = data.get('msaDatabases', '')
            msa_db.remove("on")
            msa_db_str = json.dumps(msa_db)  # ← JSON 문자열로 변환
            print(msa_db_str)
            
            subprocess.run([f"/opt/anaconda3/envs/python_11/bin/python /opt/platform/smart_bench/static/tools/msa_run.py \"{protein}\" '{msa_db_str}' \"{results_dir}\""], shell=True)
            
            # msa.json이 생길 때까지 기다림
            msa_json = f"{results_dir}/msa.json"
            timeout = 200
            elapsed = 0
            while not os.path.exists(msa_json) and elapsed < timeout:
                time.sleep(1)
                elapsed += 1

            if not os.path.exists(msa_json):
                raise FileNotFoundError("msa.json이 생성되지 않았습니다.")

            # 파일이 존재하면 열기
            with open(msa_json, 'r') as f:
                msa = json.load(f)
            
            smiles = data.get('molSeq', '')
            num_mols = int(data.get('numMols', 10))
            property_name = data.get('property_name', 'QED')
            minimize = data.get('minimize', False)
            min_similarity = float(data.get('min_similarity', 0.3))
            particles = int(data.get('particles', 20))
            iterations = int(data.get('iterations', 3))
            
            poses = int(data.get('poses', 2))

            # 함수 호출
            folded_protein = nvidiaOpenfold(protein, msa, results_dir)
            molmin_response = nvidiaMolmin(smiles, num_mols, property_name, minimize, min_similarity, particles, iterations)
            generated_ligands = '\n'.join([v['sample'] for v in molmin_response])
            print(generated_ligands)
            # generated_ligands = '\n'.join([ v['smiles']for v in molmin_response['generated']])
            diffdock_response = nvidiaDiffdock(folded_protein, generated_ligands, poses)

            # print(molmin_response)
            
            context = {
                'protein_pdb': diffdock_response['protein'],
                'ligands': [
                    {
                        'pose': i + 1,
                        'smiles': molmin_response[i]['sample'],
                        'score': molmin_response[i]['score'],
                        'sdf': diffdock_response['ligand_positions'][i][0]
                    }
                    for i in range(len(diffdock_response['ligand_positions']))
                ]
            }
            
            # with open('/opt/platform/smart_bench/apps/screening/test_nvidia.json', 'r') as f:
            #     context = json.load(f)
            
            screening_result = f"{results_dir}/screening.json"
            with open(screening_result, "w") as f:
                json.dump(context, f, indent=2)
                
            return JsonResponse({"status": "done", "screening_result" : screening_result})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)
        
def generateMoleculesSample(request):
    if request.method == "POST":    
        try:
            results_dir = f"{outputPath()}/complex_sample"
            
            data = json.loads(request.body)
            protein = data.get('proteinName', '')
            smiles = data.get('smilesName', '')
            
            screening_result = f"{results_dir}/{protein}_{smiles}/screening.json"
                
            return JsonResponse({"status": "done", "screening_result" : screening_result})  # ← 클라이언트에 저장 완료만 알림
        except Exception as e:
            print("에러 발생:", e)
            return JsonResponse({"error": str(e)}, status=500)

def showMolecules(request):
    screening_result = request.GET.get('screening_result')
    print("screening_result : ", screening_result)
    with open(screening_result, 'r') as f:
        context = json.load(f)
    
    # print("showMolecules : 세션에서 가져온 결과:", context)

    if context is None:
        return JsonResponse({"error": "세션에 저장된 결과가 없습니다."}, status=404)

    return JsonResponse(context)

def screeningInputPage(request):
    sections = ["Select Sample","Protein Seqenuce", "Chemical Generator", "Chemical Library"]
    databases = ['ChEMBL', 'ZINC', 'LNP Library']
    
    ABL1_2HZI = ["GAMDPSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQES",
                  "GAMDPSPNYDKWEMERTDITMKHKLGGGQYGEVYEGVWKKYSLTVAVKTLKEDTMEVEEFLKEAAVMKEIKHPNLVQLLGVCTREPPFYIITEFMTYGNLLDYLRECNRQEVNAVVLLYMATQISSAMEYLEKKNFIHRDLAARNCLVGENHLVKVADFGLSRLMTGDTYTAHAGAKFPIKWTAPESLAYNKFSIKSDVWAFGVLLWEIATYGMSPYPGIDLSQVYELLEKDYRMERPEGCPEKVYELMRACWQWNPSDRPSFAEIHQAFETMFQES"]
    AKT1_3CQW = ["GAMDPRVTMNEFEYLKLLGKGTFGKVILVKEKATGRYYAMKILKKEVIVAKDEVAHTLTENRVLQNSRHPFLTALKYSFQTHDRLCFVMEYANGGELFFHLSRERVFSEDRARFYGAEIVSALDYLHSEKNVVYRDLKLENLMLDKDGHIKITDFGLCKEGIKDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHEKLFELILMEEIRFPRTLGPEAKSLLSGLLKKDPKQRLGGGSEDAKEIMQHRFFAGIVWQHVYEKKLSPPFKPQVTSETDTRYFDEEFTAQMITITPPDQDDSMECVDSERRPHFPQFDYSASSTAGRPRTTSFAE",
                  "GRPRTTSFAE"]
    AKT2_3D0E = ["KVTMNDFDYLKLLGKGTFGKVILVREKATGRYYAMKILRKEVIIAKDEVAHTVTESRVLQNTRHPFLTALKYAFQTHDRLCFVMEYANGGELFFHLSRERVFTEERARFYGAEIVSALEYLHSRDVVYRDIKLENLMLDKDGHIKITDFGLCKEGISDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHERLFELILMEEIRFPRTLSPEAKSLLAGLLKKDPKQRLGGGPSDAKEVMEHRFFLSINWQDVVQKKLLPPFKPQVTSEVDTRYFDDEFTAQSITITPPDRYDSLGLLELDQRTHFPQFDYSASIR",
                  "KVTMNDFDYLKLLGKGTFGKVILVREKATGRYYAMKILRKEVIIAKDEVAHTVTESRVLQNTRHPFLTALKYAFQTHDRLCFVMEYANGGELFFHLSRERVFTEERARFYGAEIVSALEYLHSRDVVYRDIKLENLMLDKDGHIKITDFGLCKEGISDGATMKTFCGTPEYLAPEVLEDNDYGRAVDWWGLGVVMYEMMCGRLPFYNQDHERLFELILMEEIRFPRTLSPEAKSLLAGLLKKDPKQRLGGGPSDAKEVMEHRFFLSINWQDVVQKKLLPPFKPQVTSEVDTRYFDDEFTAQSITITPPDRYDSLGLLELDQRTHFPQFDYSASIR"]
    
    samples = [
        {"name": "ABL1_2HZI", "seqs": ABL1_2HZI, "smiles" : "Cn1c(=O)c(-c2c(Cl)cccc2Cl)cc2cnc(Nc3ccc(C(=O)N4CCNCC4)cc3)nc21"},
        {"name": "AKT1_3CQW", "seqs": AKT1_3CQW, "smiles" : "CCc1c[nH]c2ncnc(N3CCC(N)(CNC(=O)c4ccc(F)cc4F)C3)c12"},
        {"name": "AKT2_3D0E", "seqs": AKT2_3D0E, "smiles" : "CCn1c(-c2nonc2N)nc2c(C#CC(C)(C)O)ncc(OCCCNCC(O)CO)c21"},
    ]
    
    return render(request, 'screening/screen/input.html', {
        'sections': sections,
        'databases': databases,
        'samples_json': json.dumps(samples)
    })
    
def screeningOutputPage(request):
    return render(request, 'screening/screen/result.html')