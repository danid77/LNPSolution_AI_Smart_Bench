import subprocess
import os
import sys
import shutil
import json
import pandas as pd

from datetime import datetime

def runBoltz(input_path, output_path):
    if not input_path or not output_path:
        print("Error: Both input and output paths must be provided.")
        sys.exit(1)
    
    recycling_steps = 3
    sampling_steps = 100
    num_workers = 4
    
    # Boltz 실행 명령어 구성
    run = [
        "/opt/anaconda3/envs/boltz/bin/boltz",
        "predict",
        f"{input_path}",
        "--out_dir", f"{output_path}",
        "--cache", "/opt/git_tools/boltz/model/",
        "--output_format", "pdb",
        "--recycling_steps", str(recycling_steps),
        "--sampling_steps", str(sampling_steps),
        "--diffusion_samples", "1",
        "--num_workers", str(num_workers),
        "--use_msa_server"
    ]
     # 계산 시작 시간 기록
    start_time = datetime.now()
    print(f"Boltz 계산 중... 잠시만 기다려 주세요. (시작 시간: {start_time})")

    # Boltz 실행
    result_boltz = subprocess.run(run, capture_output=True, text=True)

    # 계산 종료 시간 기록
    end_time = datetime.now()

    if result_boltz.returncode != 0:
        print("Boltz 실행 중 오류가 발생했습니다.")
        print(f"오류 메시지: {result_boltz.stderr}")
    else:
        elapsed_time = end_time - start_time
        print(f"결과는 {output_path}에서 저장되었습니다!")
        print(f"결과 파일은 boltz_results_[input 파일 이름]/predictions/[input 파일 이름] 내에 .cif 형태로 저장되었습니다!")
        print(f"총 소요 시간: {elapsed_time}")

        # 결과 디렉토리에서 스코어 추출 및 CSV 저장
        # 입력 파일 이름에서 확장자 제거
        input_basename = os.path.splitext(os.path.basename(input_path))[0]
        predictions_folder = os.path.join(output_path, f"boltz_results_{input_basename}/predictions/{input_basename}")
        
        if not os.path.exists(predictions_folder):
            print(f"예상된 결과 폴더가 없습니다: {predictions_folder}")
        else:
            score_files = [f for f in os.listdir(predictions_folder) if f.startswith("confidence_") and f.endswith(".json")]
            
            if not score_files:
                print("스코어 파일이 발견되지 않았습니다.")
            else:
                # 주요 스코어 추출
                data = []
                for score_file in score_files:
                    score_path = os.path.join(predictions_folder, score_file)
                    with open(score_path, "r") as f:
                        json_data = json.load(f)
                    
                    # 주요 스코어 추출
                    confidence_score = json_data.get("confidence_score")
                    ptm = json_data.get("ptm")
                    iptm = json_data.get("iptm")
                    ligand_iptm = json_data.get("ligand_iptm")
                    complex_plddt = json_data.get("complex_plddt")
                    complex_iplddt = json_data.get("complex_iplddt")
                    
                    # 데이터 추가
                    data.append({
                        "File": score_file,
                        "Confidence Score": confidence_score,
                        "PTM": ptm,
                        "IPTM": iptm,
                        "Ligand IPTM": ligand_iptm,
                        "Complex pLDDT": complex_plddt,
                        "Complex ipLDDT": complex_iplddt,
                    })

                # DataFrame으로 변환 및 반올림
                df = pd.DataFrame(data)
                df = df.round(3)  # 소수점 세 자리로 반올림

                # CSV 저장
                csv_path = os.path.join(predictions_folder, "scores_summary.csv")
                df.to_csv(csv_path, index=False)
                print(f"스코어 요약이 CSV 파일로 저장되었습니다: {csv_path}")
    
    return predictions_folder


if __name__ == "__main__":    
    if len(sys.argv) < 2:
        sys.exit(1)
        
    input_fasta = sys.argv[1].strip('"').strip()
    out_path = sys.argv[2].strip('"').strip()
    
    runBoltz(input_fasta, out_path)