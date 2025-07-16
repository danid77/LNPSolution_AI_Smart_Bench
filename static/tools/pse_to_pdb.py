import os
from argparse import ArgumentParser
from pymol import cmd
import __main__

parser = ArgumentParser()
parser.add_argument("--pse", required=True)
args = parser.parse_args()
covert_pdb = args.pse.split(".")[0] + ".pdb"


__main__.pymol_argv = ['pymol', '-qc']  # quiet, no GUI
import pymol
pymol.finish_launching()

cmd.load(args.pse)           # ① .pse 파일 불러오기
cmd.save(covert_pdb, "all")          # ② 모든 구조를 .pdb로 저장
cmd.quit()

os.remove(args.pse)
