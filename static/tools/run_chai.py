import subprocess
import os
import sys
import shutil

def cleanAndListCifFiles(cc_dir):
    for f in os.listdir(cc_dir):
        file_path = os.path.join(cc_dir, f)
        # 파일인 경우
        if os.path.isfile(file_path):
            if not f.lower().endswith('.cif'):
                os.remove(file_path)
                print(f"파일 삭제: {file_path}")
        # 폴더인 경우
        elif os.path.isdir(file_path):
            shutil.rmtree(file_path)
            print(f"폴더 삭제: {file_path}")
            
    # 삭제 후, 남은 .cif 파일 리스트 작성
    cif_files = [os.path.join(cc_dir, f) for f in os.listdir(cc_dir) if f.lower().endswith('.cif')]
    return cif_files

def runChaiLabFold(input_path, output_path):
    env = os.environ.copy()
    env["MPLBACKEND"] = "Agg"
    env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"  # 메모리 조각화 방지

    try:
        process = subprocess.run(
            [
                "/opt/anaconda3/envs/chai/bin/chai-lab", "fold",
                input_path,
                output_path,
                "--device", "cuda:1",
                "--use-esm-embeddings",
                "--use-msa-server",
                "--num-diffn-samples", "1",
                "--num-diffn-timesteps", "50",
                "--low-memory"
            ],
            check=True,
            stdout=sys.stdout,
            stderr=sys.stderr,
            env=env
        )
        print("Chai-lab 실행 성공!")
        # cif_files = cleanAndListCifFiles(output_path)
        # for cif_file in cif_files:
        #     cc_pdb = "".join(cif_file.split('.')[:-1]) + ".pdb"
        #     # print(cc_pdb)
            
        #     covert = ["obabel", cif_file, "-O", cc_pdb]
        #     subprocess.run(covert)
        #     print(f"파일 삭제: {cif_file}")
        #     os.remove(cif_file)
            
        cif_file = [os.path.join(output_path, f) for f in os.listdir(output_path) if f.lower().endswith('.cif')][0]
        cc_pdb = "".join(cif_file.split('.')[:-1]) + ".pdb"
        
        covert = ["obabel", cif_file, "-O", cc_pdb]
        subprocess.run(covert)
    except subprocess.CalledProcessError as e:
        print("Chai-lab 실행 실패!")
        print(e.stderr)


if __name__ == "__main__":    
    if len(sys.argv) < 2:
        sys.exit(1)
        
    input_fasta = sys.argv[1].strip('"').strip()
    out_path = sys.argv[2].strip('"').strip()
    
    runChaiLabFold(input_fasta, out_path)