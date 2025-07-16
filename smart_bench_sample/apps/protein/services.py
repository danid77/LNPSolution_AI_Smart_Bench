from pathlib import Path
import json, os
import subprocess
import time
import requests
import shutil
import stat

from apps import getKey

def on_rm_error(func, path, exc_info):
    # 권한 바꿔서 다시 삭제 시도
    os.chmod(path, stat.S_IWRITE)
    func(path)

def setSample():
    sequences1 = "MATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQMATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ"
    sequences2 = "MLARALLLCAVLALSHTANPCCSHPCQNRGVCMSVGFDQYKCDCTRTGFYGENCSTPEFLTRIKLFLKPTPNTVHYILTHFKGFWNVVNNIPFLRNAIMSYVLTSRSHLIDSPPTYNADYGYKSWEAFSNLSYYTRALPPVPDDCPTPLGVKGKKQLPDSNEIVEKLLLRRKFIPDPQGSNMMFAFFAQHFTHQFFKTDHKRGPAFTNGLGHGVDLNHIYGETLARQRKLRLFKDGKMKYQIIDGEMYPPTVKDTQAEMIYPPQVPEHLRFAVGQEVFGLVPGLMMYATIWLREHNRVCDVLKQEHPEWGDEQLFQTSRLILIGETIKIVIEDYVQHLSGYHFKLKFDPELLFNKQFQYQNRIAAEFNTLYHWHPLLPDTFQIHDQKYNYQQFIYNNSILLEHGITQFVESFTRQIAGRVAGGRNVPPAVQKVSQASIDQSRQMKYQSFNEYRKRFMLKPYESFEELTGEKEMSAELEALYGDIDAVELYPALLVEKPRPDAIFGETMVEVGAPFSLKGLMGNVICSPAYWKPSTFGGEVGFQIINTASIQSLICNNVKGCPFTSFSVPDPELIKTVTINASSSRSGLDDINPTVLLKERSTEL"
    sequences3 = "MSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIALMSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIALMSTESMIRDVELAEEALPKKTGGPQGSRRCLFLSLFSFLIVAGATTLFCLLHFGVIGPQREEFPRDLSLISPLAQAVRSSSRTPSDKPVAHVVANPQAEGQLQWLNRRANALLANGVELRDNQLVVPSEGLYLIYSQVLFKGQGCPSTHVLLTHTISRIAVSYQTKVNLLSAIKSPCQRETPEGAEAKPWYEPIYLGGVFQLEKGDRLSAEINRPDYLDFAESGQVYFGIIAL"
    sequences4 = "MADSRDPASDQMQHWKEQRAAQKADVLTTGAGNPVGDKLNVITVGPRGPLLVQDVVFTDEMAHFDRERIPERVVHAKGAGAFGYFEVTHDITKYSKAKVFEHIGKKTPIAVRFSTVAGESGSADTVRDPRGFAVKFYTEDGNWDLVGNNTPIFFIRDPILFPSFIHSQKRNPQTHLKDPDMVWDFWSLRPESLHQVSFLFSDRGIPDGHRHMNGYGSHTFKLVNANGEAVYCKFHYKTDQGIKNLSVEDAARLSQEDPDYGIRDLFNAIATGKYPSWTFYIQVMTFNQAETFPFNPFDLTKVWPHKDYPLIPVGKLVLNRNPVNYFAEVEQIAFDPSNMPPGIEASPDKMLQGRLFAYPDTHRHRLGPNYLHIPVNCPYRARVANYQRDGPMCMQDNQGGAPNYYPNSFGAPEQQPSALEHSIQYSGEVRRFNTANDDNVTQVRAFYVNVLNEEQRKRLCENIAGHLKDAQIFIQKKAVKNFTEVHPDYGSHIQALLDKYNAEKPKNAIHTFVQSGSHLAAREKANL"
    sequences5 = "MDGETAEEQGGPVPPPVAPGGPGLGGAPGGRREPKKYAVTDDYQLSKQVLGLGVNGKVLECFHRRTGQKCALKLLYDSPKARQEVDHHWQASGGPHIVCILDVYENMHHGKRCLLIIMECMEGGELFSRIQERGDQAFTEREAAEIMRDIGTAIQFLHSHNIAHRDVKPENLLYTSKEKDAVLKLTDFGFAKETTQNALQTPCYTPYYVAPEVLGPEKYDKSCDMWSLGVIMYILLCGFPPFYSNTGQAISPGMKRRIRLGQYGFPNPEWSEVSEDAKQLIRLLLKTDPTERLTITQFMNHPWINQSMVVPQTPLHTARVLQEDKDHWDEVKEEMTSALATMRVDYDQVKIKDLKTSNNRLLNKRRKKQAGSSSASQGCNNQ"
    
    # 샘플 딕셔너리로 구성
    samples = [
        {"name": "SODC", "seq": sequences1},
        {"name": "COX2", "seq": sequences2},
        {"name": "TNFA", "seq": sequences3},
        {"name": "CATA", "seq": sequences4},
        {"name": "MAPK3", "seq": sequences5},
    ]
    return samples

def complexSample():
    sequences1 = ["PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF", "PQVTLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGATLNF"]
    sequences2 = ["DQPMEEEEVETFAFQAEIAQLMSLIINTFYSNKEIFLRELISNSSDALDKIRYETLTDPSKLDSGKELHINLIPNKQDRTLTIVDTGIGMTKADLINNLGTIAKSGTKAFMEALQAGADISMIGQFGVGFYSAYLVAEKVTVITKHNDDEQYAWESSAGGSFTVRTDTGEPMGRGTKVILHLKEDQTEYLEERRIKEIVKKHSQFIGYPITLFVEKERDKEVSDDEAE"]
    sequences3 = ["NPCCSHPCQNRGVCMSVGFDQYKCDCTRTGFYGENCSTPEFLTRIKLFLKPTPNTVHYILTHFKGFWNVVNNIPFLRNAIMSYVLTSRSHLIDSPPTYNADYGYKSWEAFSNLSYYTRALPPVPDDCPTPLGVKGKKQLPDSNEIVEKLLLRRKFIPDPQGSNMMFAFFAQHFTHQFFKTDHKRGPAFTNGLGHGVDLNHIYGETLARQRKLRLFKDGKMKYQIIDGEMYPPTVKDTQAEMIYPPQVPEHLRFAVGQEVFGLVPGLMMYATIWLREHNRVCDVLKQEHPEWGDEQLFQTSRLILIGETIKIVIEDYVQHLSGYHFKLKFDPELLFNKQFQYQNRIAAEFNTLYHWHPLLPDTFQIHDQKYNYQQFIYNNSILLEHGITQFVESFTRQIAGRVAGGRNVPPAVQKVSQASIDQSRQMKYQSFNEYRKRFMLKPYESFEELTGEKEMSAELEALYGDIDAVELYPALLVEKPRPDAIFGETMVEVGAPFSLKGLMGNVICSPAYWKPSTFGGEVGFQIINTASIQSLICNNVKGCPFTSFSVP", "NPCCSHPCQNRGVCMSVGFDQYKCDCTRTGFYGENCSTPEFLTRIKLFLKPTPNTVHYILTHFKGFWNVVNNIPFLRNAIMSYVLTSRSHLIDSPPTYNADYGYKSWEAFSNLSYYTRALPPVPDDCPTPLGVKGKKQLPDSNEIVEKLLLRRKFIPDPQGSNMMFAFFAQHFTHQFFKTDHKRGPAFTNGLGHGVDLNHIYGETLARQRKLRLFKDGKMKYQIIDGEMYPPTVKDTQAEMIYPPQVPEHLRFAVGQEVFGLVPGLMMYATIWLREHNRVCDVLKQEHPEWGDEQLFQTSRLILIGETIKIVIEDYVQHLSGYHFKLKFDPELLFNKQFQYQNRIAAEFNTLYHWHPLLPDTFQIHDQKYNYQQFIYNNSILLEHGITQFVESFTRQIAGRVAGGRNVPPAVQKVSQASIDQSRQMKYQSFNEYRKRFMLKPYESFEELTGEKEMSAELEALYGDIDAVELYPALLVEKPRPDAIFGETMVEVGAPFSLKGLMGNVICSPAYWKPSTFGGEVGFQIINTASIQSLICNNVKGCPFTSFSVP"]
    sequences4 = ["SGFRKMAFPSGKVEGCMVQVTCGTTTLNGLWLDDVVYCPRHVICTSEDMLNPNYEDLLIRKSNHNFLVQAGNVQLRVIGHSMQNCVLKLKVDTANPKTPKYKFVRIQPGQTFSVLACYNGSPSGVYQCAMRPNFTIKGSFLNGSCGSVGFNIDYDCVSFCYMHHMELPTGVHAGTDLEGNFYGPFVDRQTAQAAGTDTTITVNVLAWLYAAVINGDRWFLNRFTTTLNDFNLVAMKYNYEPLTQDHVDILGPLSAQTGIAVLDMCASLKELLQNGMNGRTILGSALLEDEFTPFDVVRQCSGVTFQ", "XAVLXX"]
    sequences5 = ["GPEFEAAACKKYMSKLRTIVAAQSRFLSTYDGAENLCLEDIYTENTLEVRTEVGMAGPLHKSPAALGLEELFSPNGHLNEDADTVLVVGEAGSGKSTLLQQVHLLWATGQDFQEFLFVFPFSCRQLQCVARPLSVMTLLFEHCCWPDVGQQDVFQFLLDHPDRILLTFDGFDEFKFKFTDHERHCSPTDPTSVQTLLFNLLQGNLLKNARKVLTSRPDAVSAFLRKYVRTEFNLKGFSEEGIELYLRKCHREPGVADRLIHLLQTTSALHGLCHLPVFSWMVSKCHQELLLQDGGSPKTTTDMYLLILQHFLRHASLPDSASQGLGPSLLQGRLPTLLRLGQLALWGLGMCCYVFSAQQLQAAQVDPDDISLGFLVQAQGVVPGSTAPLEFLHITFQCFLAAFYLVLSTDVPTASLRYLFNCRRPGSSPLSRLLPRLCVQGSEHKESTVAALLQKTEPHNLQITAAFLAGLLSREHRDLLAACQASERSLLRRRACARWCLARSLHKHFRSIPPAVPGEAKSMHAMPGFLWLIRSLYEMQEERLAQEAVRGLNVEHLKLTFCGVGPAECAALAFVLRHLRRPVALQLDHNSVGDIGVEQLLPCLGACKALYLRDNNISDRGICKLIEHALHCEQLQKLALFNNKLTDGCAHSVAQLLACKQNFLALRLGNNHITAEGAQVLAEGLRDNSSLQFLGFWGNKVGDKGAQALAEALSDHQSLKWLSLVGNNIGSVGAQALASMLEKNVALEELCLEENHLQDAGVCSLAEGLKRNSSLKVLKLSNNCITFVGAEALLQALASNDTILEVWLRGNPFSPEEMEALSHRDSRLLL"]

    # 샘플 딕셔너리로 구성
    samples = [
        {"name" : "HIV-1 protease", "pdb_id" : "1HVR", "seq": sequences1, "smiles" : "CC(C)C[C@@H](NC(=O)[C@@H](Cc1ccccc1)NC(=O)[C@@H](Cc2ccccc2)NC(=O)C(C)(C)C)C(=O)N[C@@H](Cc3ccccc3)C(=O)N(C)C"},
        {"name" : "Geldanamycin", "pdb_id" : "1YET", "seq": sequences2, "smiles" : "CC1CC(C2C(C1C(=O)NC3=CC=CC=C3)C(=O)C(=C(C2=O)OC)C)C(=O)C4=CC=CC=C4"},
        {"name" : "Cyclooxygenase-2 (COX-2)", "pdb_id" : "5IKR", "seq": sequences3, "smiles" : "CC1=C(C(=CC=C1)NC2=CC=CC=C2C(=O)O)C"},
        {"name" : "Covid 19", "pdb_id" : "6LU7", "seq": sequences4, "smiles" : "CC(C)C[C@H](NC(=O)[C@H](CC1=CC=CC=C1)NC(=O)[C@H](CC2=CC=CC=C2)NC(=O)C(C)(C)C)C(=O)N[C@H](CC3=CC=CC=C3)C(=O)N(C)C"},
        {"name" : "NOD2 gene of Crohn's disease", "pdb_id" : "5IRN", "seq": sequences5, "smiles" : "c1nc(c2c(n1)n(cn2)[C@H]3[C@@H]([C@@H]([C@H](O3)CO[P@](=O)(O)OP(=O)(O)O)O)O)N"},
    ]
    return samples


def proteinMakeResultFolder(main_results_dir):
    openfold_result_dir = f"{main_results_dir}/openfold"
    alphafold_result_dir = f"{main_results_dir}/alphafold"
    esmfold_result_dir = f"{main_results_dir}/esmfold"
    boltz_result_dir = f"{main_results_dir}/boltz"
    chai_result_dir = f"{main_results_dir}/chai"

    for dir in [openfold_result_dir, alphafold_result_dir, esmfold_result_dir, boltz_result_dir, chai_result_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir, exist_ok=True)
    
    return openfold_result_dir, alphafold_result_dir, esmfold_result_dir, boltz_result_dir, chai_result_dir

def runDownloadPDB(results_dir, sample_name):
    print("input sample : ", sample_name)
    
    sample_path = "/opt/platform/smart_bench/static/output/protein_sample/origin_pdb"
    
    # 원본 파일 경로
    source_file = os.path.join(sample_path, f"{sample_name}_origin.pdb")
    destination_file = os.path.join(results_dir, "origin.pdb")

    # 복사
    if os.path.isfile(source_file):
        shutil.copy2(source_file, destination_file)  # 메타데이터까지 유지하려면 copy2 사용
    else:
        print("Response : fail", "500")

def runOpenfold(results_dir, msa_db, protein):
    msa_db_str = json.dumps(msa_db)  # ← JSON 문자열로 변환
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
        
    key = getKey()
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

    # ---------------------------------------------------------
    # View response
    # ---------------------------------------------------------
    if response.status_code == 200:
        output_file.write_text(response.text)
        print(f"Response output to file: {output_file}")
        return "############## OpenFold2 Complete ###############"
    else:
        print(f"Unexpected HTTP status: {response.status_code}")
        print(f"Response: {response.text}") 
        return "############## OpenFold2 Fail ###############"


# def runAlphafold2(results_dir, protein):
#     url = os.getenv("URL", "https://health.api.nvidia.com/v1/biology/deepmind/alphafold2")
#     status_url = os.getenv("STATUS_URL", "https://health.api.nvidia.com/v1/status")

#     sequence = (protein)
#     output_file = Path(f"{results_dir}/af2_output.json")

#     # Initial request
#     headers = {
#         "content-type": "application/json",
#         "Authorization": f"Bearer {getKey()}",
#         "NVCF-POLL-SECONDS": "300",
#     }
#     data = {
#         "sequence": sequence,
#         "algorithm": "mmseqs2",
#         "e_value": 0.0001,
#         "iterations": 1,
#         "databases": ["small_bfd"],
#         "relax_prediction": False,
#         "skip_template_search" : True
#     }

#     print("Making request...")
#     response = requests.post(url, headers=headers, json=data)

#     # Check the status code
#     if response.status_code == 200:
#         output_file.write_text(response.text)
#         print(f"Response output to file: {output_file}")
#     elif response.status_code == 202:
#         print("Request accepted...")
#         # Extract reqId header
#         req_id = response.headers.get("nvcf-reqid")

#         # Poll the /status endpoint
#         while True:
#             print("Polling for response...")
#             status_response = requests.get(f"{status_url}/{req_id}", headers=headers)

#             if status_response.status_code != 202:
#                 output_file.write_text(status_response.text)
#                 print(f"Response output to file: {output_file}")
#                 break
#     else:
#         print(f"Unexpected HTTP status: {response.status_code}")
#         print(f"Response: {response.text}")

def runAlphafold2(results_dir, protein, sample_name):
    print("input sequence : ", protein)
    print("Making request...")
    
    sample_path = "/opt/platform/smart_bench/static/output/protein_sample/af2_result"
    
    # 원본 파일 경로
    source_file = os.path.join(sample_path, f"{sample_name}.json")
    destination_file = os.path.join(results_dir, f"{sample_name}.json")
    time.sleep(10)
    print("Request accepted...")
    for i in range(10):
        print("Polling for response...")
        time.sleep(0.5)
    # 복사
    if os.path.isfile(source_file):
        shutil.copy2(source_file, destination_file)  # 메타데이터까지 유지하려면 copy2 사용
        print(f"Response output to file : {results_dir}/{sample_name}.json")
    else:
        print("Response : fail", "500")

    return "############## runAlphafold2 Complete ##############"

def runEsmfold(results_dir, protein):
    invoke_url = "https://health.api.nvidia.com/v1/biology/nvidia/esmfold"

    headers = {
        "Authorization": f"Bearer {getKey()}",
        "Accept": "application/json",
    }

    payload = {
    "sequence": protein
    }

    # re-use connections
    session = requests.Session()
    print(session)
    response = session.post(invoke_url, headers=headers, json=payload)
    print(response)
    response.raise_for_status()
    
    if response.status_code == 200:
        response_body = response.json()
        print(response_body)
        with open(f"{results_dir}/esmfold.json", "w") as f:
            json.dump(response_body, f, indent=2)

    else:
        print(f"Unexpected HTTP status: {response.status_code}")
        print(f"Response: {response.text}")
    return "############## runEsmfold Complete ##############"

# boltz 인풋파일 생성 (단백질 서열 1개만 처리, 리간드는 제외)
def runBoltzChai(results_dir, boltz_result_dir, chai_result_dir, protein):
    boltz_input_file = f"{results_dir}/boltz_input.fasta"
    
    print("boltz_input_file")
    print(protein)

    with open(boltz_input_file, "w") as fasta:
        fasta.write(f">A|protein|\n")
        fasta.write(f"{protein}\n")
    
    subprocess.run([f"CUDA_VISIBLE_DEVICES=1 /opt/anaconda3/envs/boltz/bin/python /opt/platform/smart_bench/static/tools/run_boltz.py {boltz_input_file} {boltz_result_dir}"], shell=True)
    print("############## Boltz Complete ##############")
    
    time.sleep(10)
    
    chai_input_file = f"{results_dir}/chai_input.fasta"
    
    print("chai_input_file")
    print(protein)

    with open(chai_input_file, "w") as fasta:
        fasta.write(f">protein|name=chain_A\n")
        fasta.write(f"{protein}\n")
    
    subprocess.run([f"/opt/anaconda3/envs/chai/bin/python /opt/platform/smart_bench/static/tools/run_chai.py {chai_input_file} {chai_result_dir}"], shell=True)
    print("############## Chai Complete ##############")
    
    os.remove(boltz_input_file)
    os.remove(chai_input_file)
    
    return "############## runBoltzChai Complete ##############"
   
# def runBoltz(input_fasta, results_dir):

# # chai 인풋파일 생성 (단백질 서열 1개만 처리, 리간드는 제외)
# def runChai(results_dir, chai_result_dir, protein):
#     chai_input_file = f"{results_dir}/chai_input.fasta"
    
#     print("chai_input_file")
#     print(protein)

#     with open(chai_input_file, "w") as fasta:
#         fasta.write(f">protein|name=chain_A\n")
#         fasta.write(f"{protein}\n")
    
#     subprocess.run([f"/opt/anaconda3/envs/chai/bin/python /opt/platform/smart_bench/static/tools/run_chai.py {chai_input_file} {chai_result_dir}"], shell=True)

def proteinAlignment(results_dir):
    origin_pdb = os.path.join(results_dir, "origin.pdb")

    tools = {
        "openfold": ("openfold/openfold.json", "json", lambda d: d["structures_in_ranked_order"][0]["structure"]),
        "alphafold": (lambda d: f"alphafold/{os.listdir(os.path.join(d, 'alphafold'))[0]}", "json", lambda d: d[0]),
        "esmfold": ("esmfold/esmfold.json", "json", lambda d: d["pdbs"][0]),
        "boltz": ("boltz/boltz_results_boltz_input/predictions/boltz_input/boltz_input_model_0.pdb", "text", None),
        "chai": ("chai/predmodel_idx_0.pdb", "text", None),
    }

    pdb_paths = {}

    for tool, (path_def, file_type, extractor) in tools.items():
        src_path = os.path.join(results_dir, path_def(results_dir) if callable(path_def) else path_def)
        out_path = os.path.join(results_dir, f"{tool}.pdb")
        pdb_paths[tool] = out_path

        with open(src_path, "r") as src_file, open(out_path, "w") as out_file:
            if file_type == "json":
                json_data = json.load(src_file)
                out_file.write(extractor(json_data))
            else:
                # text 파일은 스트리밍 복사
                for line in src_file:
                    out_file.write(line)

        # 결과 디렉터리 삭제
        dir_path = os.path.join(results_dir, tool)
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path, onerror=on_rm_error)
            print(f"폴더 삭제: {dir_path}")
        else:
            print(f"폴더 없음: {dir_path}")

    # 구조 정렬 실행
    align_cmd = ["/opt/schrodinger2024-2/utilities/structalign", origin_pdb] + list(pdb_paths.values())
    align = subprocess.run(align_cmd, capture_output=True, text=True, cwd=results_dir)

    # 임시 PDB 파일 삭제
    os.remove(origin_pdb)
    for path in pdb_paths.values():
        if os.path.exists(path):
            os.remove(path)
            print(f"파일 삭제: {path}")
            
    return align.stdout