import tempfile
import json
import os
import subprocess
from datetime import datetime
import pandas as pd

def molsparkProcess(smiles, model_type, strategy, temperature, num_samples):
    # print(smiles, model_type, strategy, temperature, num_samples)
    current_time = datetime.now().strftime("%Y%m%d_%Hh_%Mm_%Ss")
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
            "model_file": f"/opt/git-tools/REINVENT4/models/{model_type}.prior",
            "smiles_file": smiles_path,
            "sample_strategy": strategy,
            "output_file": output_file,
            "num_smiles": int(num_samples),
            "unique_molecules": True,
            "randomize_smiles": True
        }
    }
    
    # temperature가 None이 아니면 추가
    if temperature is not None:
        reinvent_option["parameters"]["temperature"] = temperature

    # 임시 JSON 파일 생성
    with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".json") as option_json:
        json.dump(reinvent_option, option_json, indent=4)
        option_json_path = option_json.name

    print(f"임시 JSON 파일 경로: {option_json_path}")
    print(output_file)
    
    cal_output_dir = f"/opt/platform/smart_bench/static/output/molspark/{current_time}"
    if not os.path.exists(cal_output_dir):
        os.mkdir(cal_output_dir)
    
    reinvent = ["/opt/anaconda3/envs/reinvent4/bin/reinvent", "-f", "json", option_json_path]
    reinvent_cal = ["/opt/anaconda3/envs/reinvent4/bin/python", 
                    "/opt/git-tools/REINVENT4/reinvent/notebooks/property_cal2.py", 
                    output_file, cal_output_dir]
    process1 = subprocess.run(reinvent, capture_output=True, text=True)
    process2 = subprocess.run(reinvent_cal, capture_output=True, text=True)
    print("STDOUT:", process1.stdout)
    print("STDERR:", process1.stderr)
    print("STDOUT:", process2.stdout)
    print("STDERR:", process2.stderr)
    
    os.remove(smiles_path)
    os.remove(option_json_path)
    os.remove(output_file)
    
    csv_list = os.listdir(cal_output_dir)
    csv_file = [file for file in csv_list if file.find(".csv") !=-1][0]
    cal_output_file = f"{cal_output_dir}/{csv_file}"
    
    return cal_output_file

def csv_to_dict(_csv_file):
    """ CSV 파일을 Pandas DataFrame으로 읽어 딕셔너리 리스트로 변환 """
    try:
        csv_file = pd.read_csv(_csv_file).astype(object)
        csv_file = csv_file.where(pd.notnull(csv_file), None)  # 🔥 `_csv_file` 대신 `csv_file` 사용!
        return csv_file.to_dict(orient='records')
    except Exception as e:
        print(f"CSV 변환 오류: {e}")
        return None