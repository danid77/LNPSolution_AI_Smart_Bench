import subprocess
from argparse import ArgumentParser
import os
import re
import pandas as pd

parser = ArgumentParser()
parser.add_argument("-r", "--recepter", required=True)
parser.add_argument("-l", "--ligand", required=True)
# parser.add_argument("-c", "--CNNmodel", default="crossdock_default2018")
# parser.add_argument("-o", "--output", default=f"/home/dwpdp/test/gnina/results/.sdf")
args = parser.parse_args()

data = []
lines_list = []

def name_export(path):
    base_name = os.path.basename(path)
    file_name, _ = os.path.splitext(base_name)
    return file_name

command = f"/opt/platform/smart_bench/static/tools/gnina \
            -r {args.recepter} \
            -l {args.ligand} \
            --score_only \
            --log {name_export(args.ligand)}.txt \
            --device 2"

print(f"Running: {command}")
gnina_score = subprocess.run(command, shell=True, capture_output=True)
# # print(gnina_score.stdout)

gnina_score_output = str(gnina_score.stdout)

log_entries = gnina_score_output.split("##")
print(log_entries[1:3])

for n in range(len(log_entries[1:])):
    ligand_match = log_entries[n + 1].split()[0] if log_entries[n + 1].split()[0] else None
    affinity_match = re.search(r"Affinity: (-?\d+\.\d+)", log_entries[n])
    cnn_score_match = re.search(r"CNNscore: (\d+\.\d+)", log_entries[n])
    cnn_affinity_match = re.search(r"CNNaffinity: (-?\d+\.\d+)", log_entries[n])
    intra_energy_match = re.search(r"Intramolecular energy: (-?\d+\.\d+)", log_entries[n])

    if ligand_match and affinity_match and cnn_score_match and cnn_affinity_match and intra_energy_match:
        ligand = ligand_match
        affinity = float(affinity_match.group(1))
        cnn_score = float(cnn_score_match.group(1))
        cnn_affinity = float(cnn_affinity_match.group(1))
        intra_energy = float(intra_energy_match.group(1))
    
        data.append([ligand, affinity, cnn_score, cnn_affinity, intra_energy])

file_path = os.getcwd() + f"/{name_export(args.ligand)}.txt"

os.remove(file_path)
# DataFrame으로 변환
df = pd.DataFrame(data, columns=["Ligand", "Affinity_(kcal/mol)", "CNNscore", "CNNaffinity", "Intramolecular_energy"])

# DataFrame을 CSV로 저장
df.to_csv(f"{name_export(args.recepter)}_single.csv", index=False)

# print(f"Results saved to {output_file}")
