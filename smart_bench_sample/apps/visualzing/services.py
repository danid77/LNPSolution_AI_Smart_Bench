from datetime import datetime
import os, re, json
import subprocess, shutil
import pandas as pd
import xml.etree.ElementTree as ET

def visualzingFolderGenerator(main_results_dir):
    plip_result_dir = f"{main_results_dir}/plip"
    img2d_result_dir = f"{main_results_dir}/img2d"
    gnina_result_dir = f"{main_results_dir}/gnina"

    for dir in [plip_result_dir, img2d_result_dir, gnina_result_dir]:
        os.makedirs(dir, exist_ok=True)
    
    return plip_result_dir, img2d_result_dir, gnina_result_dir

def runPlip(main_file_path, plip_result_dir):
    plip_file_path = f"{plip_result_dir}/{os.path.basename(main_file_path)}"
    shutil.copy2(main_file_path, plip_file_path)
    
    # Docker로 PLIP 실행
    try:
        run_plip = subprocess.run(
            ["docker", "run", "--rm",
            "-v", f"{plip_result_dir}:/results",
            "-w", "/results",
            "pharmai/plip:latest", "-f", f"/results/{os.path.basename(main_file_path)}", "-pxty"],
            capture_output=True, text=True, check=True
        )
        
        # XML 파일 파싱
        tree = ET.parse(f'{plip_result_dir}/report.xml')
        root = tree.getroot()

        # Ligand Identifier 정보 추출
        identifiers = root.find(".//bindingsite/identifiers")
        ligand_info = {
            "longname": identifiers.findtext("longname"),
            "ligtype": identifiers.findtext("ligtype"),
            "hetid": identifiers.findtext("hetid"),
            "chain": identifiers.findtext("chain"),
            "position": identifiers.findtext("position"),
            "composite": identifiers.findtext("composite"),
            "smiles": identifiers.findtext("smiles"),
            "inchikey": identifiers.findtext("inchikey").strip(),
            "members": [m.text for m in identifiers.findall(".//member")]
        }

        # 수소결합 정보 추출
        hbonds = []
        for hbond in root.findall(".//hydrogen_bond"):
            data = {
                'interactions':"hydrogen_bond",
                'resnr': hbond.findtext('resnr'),
                'restype': hbond.findtext('restype'),
                'reschain': hbond.findtext('reschain'),
                'resnr_lig': hbond.findtext('resnr_lig'),
                'restype_lig': hbond.findtext('restype_lig'),
                'reschain_lig': hbond.findtext('reschain_lig'),
                'sidechain': hbond.findtext('sidechain'),
                'dist_h-a': hbond.findtext('dist_h-a'),
                'dist_d-a': hbond.findtext('dist_d-a'),
                'don_angle': hbond.findtext('don_angle'),
                'protisdon': hbond.findtext('protisdon'),
                'donoridx': hbond.findtext('donoridx'),
                'donortype': hbond.findtext('donortype'),
                'acceptoridx': hbond.findtext('acceptoridx'),
                'acceptortype': hbond.findtext('acceptortype'),
            }
            hbonds.append(data)

        # π-stacking 정보 추출
        pi_stacks = []
        for stack in root.findall(".//pi_stack"):
            pi_data = {
                'interactions':"pi_stack",
                "resnr": stack.findtext("resnr"),
                "restype": stack.findtext("restype"),
                "reschain": stack.findtext("reschain"),
                "resnr_lig": stack.findtext("resnr_lig"),
                "restype_lig": stack.findtext("restype_lig"),
                "reschain_lig": stack.findtext("reschain_lig"),
                "centdist": stack.findtext("centdist"),
                "angle": stack.findtext("angle"),
                "offset": stack.findtext("offset"),
                "type": stack.findtext("type"),
            }
            pi_stacks.append(pi_data)

        # DataFrame으로 변환
        print(ligand_info)
        hy_df = pd.DataFrame(hbonds)
        pi_df = pd.DataFrame(pi_stacks)
        plip_df = pd.concat([hy_df, pi_df])
        plip_df.to_csv(f"{plip_result_dir}/plip.csv", index=False)
        print(plip_df)
        
        with open(f"{plip_result_dir}/ligand_info.json", "w") as f:
            json.dump(ligand_info, f, indent=2)
        
        
        # pdb_files = [
        #     f for f in os.listdir(plip_result_dir)
        #     if f.lower().endswith('.pdb')
        # ]
    
        # pse_files = [
        #     f for f in os.listdir(plip_result_dir)
        #     if f.lower().endswith('.pse')
        # ]
        # for ipdb in range(len(pdb_files)):
        #     os.remove(f"{plip_result_dir}/{pdb_files[ipdb]}")
        
        # for ipse in range(len(pse_files)):
        #     pse_to_pdb = subprocess.run(["/opt/anaconda3/envs/pymol_env/bin/python", "/opt/platform/smart_bench/static/tools/pse_to_pdb.py", 
        #                                  "--pse", f"{plip_result_dir}/{pse_files[ipse]}"], 
        #                                 capture_output=True, text=True, check=True)
            
    except subprocess.CalledProcessError as e:
        print({
            'error': 'PLIP 분석 중 오류 발생',
            'details': e.stderr
        })
    # os.remove(plip_file_path)
    return print("200 : finish plip")

def run2DImg(main_file_path, main_file_name, img2d_result_dir):
    # file_name, _ = os.path.splitext(os.path.basename(main_file_path))
    maegz_file_path = f"{img2d_result_dir}/{main_file_name}.maegz"

    try:
        pdbconvert = subprocess.run(
                    ["/opt/schrodinger2024-2/utilities/pdbconvert",
                    "-ipdb", main_file_path,
                    "-omae", maegz_file_path],
                    capture_output=True, text=True, check=True
                )
    except subprocess.CalledProcessError as e:
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        raise

    try:
        os.environ["QT_QPA_PLATFORM"] = "offscreen"  # 핵심 포인트!
        ligand_interaction_diagram = subprocess.run(
            ["/opt/schrodinger2024-2/utilities/ligand_interaction_diagram",
            "-i", maegz_file_path,
            "-ligandASL", "ligand",
            "-o", f"{img2d_result_dir}/{main_file_name}_2D.png"],
            capture_output=True, text=True, check=True
        )
        os.remove(maegz_file_path)
    except subprocess.CalledProcessError as e:
        print("STDOUT:", e.stdout)
        print("STDERR:", e.stderr)
        raise
    
    return print("200 : finish 2D img")

def runGninaSocre(main_file_path, main_file_name, gnina_result_dir):
    no_ligand_pdb = f"{gnina_result_dir}/{main_file_name}_no_ligand.pdb"
    ligand_pdb = f"{gnina_result_dir}/{main_file_name}_ligand.pdb"
    ligand_sdf = f"{gnina_result_dir}/{main_file_name}_ligand.sdf"
    
    # PDB 파일에서 HETATM 정보만 추출
    with open(main_file_path, "r") as pdb:
        lines = pdb.readlines()
        ligand_lines = [line for line in lines if line.startswith("HETATM")]
        non_ligand_lines = [line for line in lines if not line.startswith("HETATM")]

    # ligand_pdb 파일 저장 (HETATM만 포함)
    with open(ligand_pdb, "w") as ligand_file:
        ligand_file.writelines(ligand_lines)

    # ligand 제외한 pdb 파일 저장 (ATOM만 포함된 pdb)
    with open(no_ligand_pdb, "w") as non_ligand_file:
        non_ligand_file.writelines(non_ligand_lines)

    # SDF 변환
    # sdf_name = ligand_pdb.split(".")[0] + ".sdf"
    export_ligand = ["/opt/openbabel/bin/obabel", "-ipdb", ligand_pdb, "-osdf", "-O", ligand_sdf]
    subprocess.run(export_ligand)
    os.remove(ligand_pdb)
    
    command = (
        f'docker run --rm --gpus "device=1" '
        f'-v {gnina_result_dir}:{gnina_result_dir} '  # 컨테이너 안에도 똑같이 마운트
        f'gnina/gnina gnina '
        f'-r {no_ligand_pdb} '
        f'-l {ligand_sdf} '
        '--score_only '
        '--device 0'
    )

    # print(f"Running: {command}")
    gnina_score = subprocess.run(command, shell=True, capture_output=True, text=True)
    # print(gnina_score.stderr)
    # print(gnina_score_output)

    gnina_score_output = gnina_score.stdout
    log_entries = gnina_score_output.split("##")
    data = []
    lines_list = []
    for n in range(len(log_entries[1:])):
        ligand_match = log_entries[n + 1].split()[0] if log_entries[n + 1].split()[0] else None
        affinity_match = re.search(r"Affinity: (-?\d+\.\d+)", log_entries[n])
        cnn_score_match = re.search(r"CNNscore: (\d+\.\d+)", log_entries[n])
        cnn_affinity_match = re.search(r"CNNaffinity: (-?\d+\.\d+)", log_entries[n])
        intra_energy_match = re.search(r"Intramolecular energy: (-?\d+\.\d+)", log_entries[n])

        if ligand_match and affinity_match and cnn_score_match and cnn_affinity_match and intra_energy_match:
            # ligand = ligand_match
            affinity = float(affinity_match.group(1))
            cnn_score = float(cnn_score_match.group(1))
            cnn_affinity = float(cnn_affinity_match.group(1))
            intra_energy = float(intra_energy_match.group(1))
        
            # data.append([ligand, affinity, cnn_score, cnn_affinity, intra_energy])
            data.append([affinity, cnn_score, cnn_affinity, intra_energy])
            
    df = pd.DataFrame(data, columns=["Affinity_(kcal/mol)", "CNNscore", "CNNaffinity", "Intramolecular_energy"])
    csv_path = f"{gnina_result_dir}/gnina_score.csv"
    df.to_csv(csv_path, index=False)
    os.remove(no_ligand_pdb)
    os.remove(ligand_sdf)
    return print("200 : finish Gnina")