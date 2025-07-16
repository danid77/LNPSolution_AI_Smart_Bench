import warnings
import os, sys
from datetime import datetime
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, AllChem, EState
import math
import umap
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.manifold import TSNE
import pandas as pd
import logging
from abc import ABC
from enum import Enum

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

logger = logging.getLogger(__name__)
warnings.simplefilter("ignore")

def calc_morgan_fingerprints(dataframe, smiles_col):
    """Calculate Morgan fingerprints on SMILES strings

    Args:
        dataframe (pd.DataFrame): dataframe containing a SMILES column for calculation

    Returns:
        pd.DataFrame: new dataframe containing fingerprints
    """
    mf = MorganFingerprint()
    fp = mf.transform(dataframe, col_name=smiles_col)  # numpy 2D array
    fp = pd.DataFrame(fp, columns=[f"fp_{i+1}" for i in range(fp.shape[1])])  # 컬럼 이름 지정
    fp.index = dataframe.index  # 원래 인덱스 유지
    return fp

class TransformationDefaults(Enum):
    MorganFingerprint = {'radius': 3, 'nBits': 2048}
    Embeddings = {}

class BaseTransformation(ABC):
    def __init__(self, **kwargs):
        self.name = None
        self.kwargs = None
        self.func = None

    def transform(self, data):
        return NotImplemented

    def transform_many(self, data):
        return list(map(self.transform, data))

    def __len__(self):
        return NotImplemented


class MorganFingerprint(BaseTransformation):
    def __init__(self, **kwargs):
        self.name = __class__.__name__.split('.')[-1]
        self.kwargs = TransformationDefaults[self.name].value
        self.kwargs.update({'radius': 3})  # ECFP6으로 변경
        self.kwargs.update(kwargs)
        self.func = AllChem.GetMorganFingerprintAsBitVect

    def transform(self, data, col_name):
        data = data[col_name]
        fp_array = []
        fail_count = 0  # 변환 실패 카운트 추가

        for mol in data:
            m = Chem.MolFromSmiles(mol)
            if m is None:
                print(f"❌ 변환 실패: {mol}")  # 문제가 있는 SMILES 확인
                fp_array.append([np.nan] * self.kwargs['nBits'])  # NaN으로 채우기
                fail_count += 1  # 실패 개수 증가
                continue

            fp = self.func(m, **self.kwargs)
            fp_array.append(np.array(fp))  # NumPy 배열로 변환
        
        print(f"\n🔹 ECFP6 : 총 {fail_count}개의 SMILES 변환 실패")  # 최종 결과 출력
        return np.vstack(fp_array)  # 모든 배열을 스택으로 정렬


def get_aromatic_proportion(mol):
    aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    return aromatic_atoms / mol.GetNumAtoms() if mol.GetNumAtoms() > 0 else 0


def numBridgeheadsAndSpiro(mol, ri=None):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro

def calculateScore(m):
    # Simplified scoring calculation
    fp = rdMolDescriptors.GetMorganFingerprint(m, 2)
    fps = fp.GetNonzeroElements()
    score1 = 0.
    nf = 0
    for bitId, v in fps.items():
        nf += v
        sfp = bitId
        score1 += v  # simplified scoring

    if nf > 0:
        score1 /= nf

    # features score
    nAtoms = m.GetNumAtoms()
    nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
    ri = m.GetRingInfo()
    nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m, ri)
    nMacrocycles = sum(1 for x in ri.AtomRings() if len(x) > 8)

    sizePenalty = nAtoms**1.005 - nAtoms
    stereoPenalty = math.log10(nChiralCenters + 1) if nChiralCenters > 0 else 0
    spiroPenalty = math.log10(nSpiro + 1) if nSpiro > 0 else 0
    bridgePenalty = math.log10(nBridgeheads + 1) if nBridgeheads > 0 else 0
    macrocyclePenalty = math.log10(2) if nMacrocycles > 0 else 0.

    score2 = - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty

    score3 = math.log(float(nAtoms) / (nf + 1)) * .5 if nf > 0 else 0.

    sascore = score1 + score2 + score3

    min = -4.0
    max = 2.5
    sascore = 11. - (sascore - min + 1) / (max - min) * 9.
    
    if sascore > 8.:
        sascore = 8. + math.log(max(sascore + 1. - 9., 1e-3))  # log 안에 0 들어가는 문제 해결
    if sascore > 10.:
        sascore = 10.0
    elif sascore < 1.:
        sascore = 1.0

    return sascore

def calculate_properties(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None,) * 14

    mol_weight = round(Descriptors.MolWt(mol), 2)
    logp = round(Descriptors.MolLogP(mol), 2)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    rot_bonds = Descriptors.NumRotatableBonds(mol)
    tpsa = round(Descriptors.TPSA(mol), 2)
    num_atoms = mol.GetNumAtoms()
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    num_aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    surface_area = round(rdMolDescriptors.CalcLabuteASA(mol), 2)
    aromatic_proportion = get_aromatic_proportion(mol)

    # Delaney LogS 예측
    logs = round(
        -0.261 * logp +
        -0.006 * mol_weight +
        -0.74 * aromatic_proportion +
        -0.006 * rot_bonds +
        0.5, 4
    )

    qed = round(QED.qed(mol), 4)
    sa_score = round(calculateScore(mol), 4)
    alogp_count = len(str(logp))

    return (
        mol_weight,         # Molecular_Weight
        logp,               # ALogP
        hbd,                # Num_H_Donors
        hba,                # Num_H_Acceptors
        rot_bonds,          # Num_RotatableBonds
        tpsa,               # Molecular_PolarSurfaceArea
        num_atoms,          # Num_Atoms
        num_rings,          # Num_Rings
        num_aromatic_rings, # Num_AromaticRings
        surface_area,       # Molecular_SurfaceArea
        logs,               # Molecular_Solubility
        qed,                # QED
        sa_score,           # SA_Score
        alogp_count,        # AlogP_Count
    )

def check_lipinski(mw, hbd, hba, logp):
    return mw <= 500.0 and hbd <= 5 and hba <= 10 and logp <= 5.0

def check_veber(tpsa, rot_bonds):
    return tpsa <= 140.0 and rot_bonds <= 10

def check_lipinski_scores(mw, hbd, hba, logp):
    """Lipinski 규칙 4가지 중 만족하는 개수를 반환"""
    score = 0
    score += mw <= 500.0
    score += hbd <= 5
    score += hba <= 10
    score += logp <= 5.0
    return score

def check_veber_scores(tpsa, rot_bonds):
    """Veber 규칙 2가지 중 만족하는 개수를 반환"""
    score = 0
    score += tpsa <= 140.0
    score += rot_bonds <= 10
    return score

def process_smiles_file(smiles_data, selected_column):
    
    # Calculate properties
    properties = smiles_data[selected_column].apply(calculate_properties)
    columns = [
        "Molecular_Weight", "ALogP", "Num_H_Donors", "Num_H_Acceptors",
        "Num_RotatableBonds", "Molecular_PolarSurfaceArea",
        "Num_Atoms", "Num_Rings", "Num_AromaticRings",
        "Molecular_SurfaceArea", "Molecular_Solubility",
        "QED", "SA_Score", "AlogP_Count"]

    properties_df = pd.DataFrame(properties.tolist(), columns=columns)

    # 규칙 확인
    properties_df['Lipinski'] = properties_df.apply(
        lambda row: check_lipinski_scores(
            row['Molecular_Weight'],
            row['Num_H_Donors'],
            row['Num_H_Acceptors'],
            row['ALogP']
        ), axis=1
    )

    properties_df['Veber'] = properties_df.apply(
        lambda row: check_veber_scores(
            row['Molecular_PolarSurfaceArea'],
            row['Num_RotatableBonds']
        ), axis=1
    )
    
    # Combine results and save
    result_df = pd.concat([smiles_data, properties_df], axis=1)
    
    # display(result_df)
    return result_df

def process_smiles_file_ecfp6(smiles_data, selected_column):
    
    # Calculate properties
    properties = smiles_data[selected_column].apply(calculate_properties)
    columns = [
        "Molecular_Weight", "ALogP", "Num_H_Donors", "Num_H_Acceptors",
        "Num_RotatableBonds", "Molecular_PolarSurfaceArea",
        "Num_Atoms", "Num_Rings", "Num_AromaticRings",
        "Molecular_SurfaceArea", "Molecular_Solubility",
        "QED", "SA_Score", "AlogP_Count"]

    properties_df = pd.DataFrame(properties.tolist(), columns=columns)

    # 규칙 확인
    properties_df['Lipinski'] = properties_df.apply(
        lambda row: check_lipinski_scores(
            row['Molecular_Weight'],
            row['Num_H_Donors'],
            row['Num_H_Acceptors'],
            row['ALogP']
        ), axis=1
    )

    properties_df['Veber'] = properties_df.apply(
        lambda row: check_veber_scores(
            row['Molecular_PolarSurfaceArea'],
            row['Num_RotatableBonds']
        ), axis=1
    )
    
    # ECFP6 계산
    fp_df = calc_morgan_fingerprints(smiles_data, selected_column)

    # Combine results and save
    result_df = pd.concat([smiles_data, properties_df, fp_df], axis=1)
    
    return result_df

def umap_coordinate(df, selected_column, cnt):
    # 2. 데이터 스케일링 (UMAP은 거리 기반이라 표준화 필요)
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df.drop(columns=[selected_column]).replace([np.inf, -np.inf], np.nan).fillna(0))

    # 다시 NaN 체크 (이 단계에서도 문제가 있으면 데이터셋 확인 필요)
    if np.isnan(X_scaled).any() or np.isinf(X_scaled).any():
        raise ValueError("X_scaled contains NaN or Inf values after preprocessing!")
    
    # 3. UMAP 적용 (1D로 차원 축소 -> UMAP1만 사용)
    reducer = umap.UMAP(n_components=int(cnt), n_neighbors=30, n_epochs=1000, metric='euclidean', random_state=42)
    X_umap = reducer.fit_transform(X_scaled)  # 결과는 (n_samples, 1) 형태
    umap_df = pd.DataFrame(X_umap, columns=[f"UMAP{i+1}" for i in range(X_umap.shape[1])])
    return umap_df

def tsne_coordinate(df, selected_column, cnt):
    """Plot t-SNE in 3D with preprocessing."""

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df.drop(columns=[selected_column]).replace([np.inf, -np.inf], np.nan).fillna(0))
    
    # 다시 NaN 체크 (이 단계에서도 문제가 있으면 데이터셋 확인 필요)
    if np.isnan(X_scaled).any() or np.isinf(X_scaled).any():
        raise ValueError("X_scaled contains NaN or Inf values after preprocessing!")

    # t-SNE 3차원 변환
    n_samples = X_scaled.shape[0]
    perplexity = min(30, (n_samples - 1) // 3)  # 안전하게 1/3 정도로 제한
    tsne = TSNE(n_components=int(cnt), perplexity=perplexity, n_iter=1000, random_state=42)
    tsne_cal = tsne.fit_transform(X_scaled)
    tsne_df = pd.DataFrame(tsne_cal, columns=[f"tSNE{i+1}" for i in range(tsne_cal.shape[1])])
    return tsne_df


if __name__ == "__main__":    
    current_time = datetime.now().strftime("%Y%m%d_%Hh%Mm%Ss")
    
    if len(sys.argv) < 2:
        print("Usage: <csv files 이름> <출력할 umap 갯수>")
        sys.exit(1)
        
    in_csv_file = sys.argv[1].strip('"').strip()
    out_csv_file_path = sys.argv[2].strip('"').strip()
    cnt = 2
    
    current_dir = os.getcwd()
    csv_file = os.path.join(current_dir, in_csv_file)
    csv_file_name = in_csv_file.split(".")[0]
    
    df = pd.read_csv(csv_file)

    # 선택할 컬럼 리스트 (우선순위: SMILES > Smiles > smiles)
    column_options = ["SMILES", "Smiles", "smiles"]

    # 존재하는 컬럼 중 하나 선택
    selected_column = next((col for col in column_options if col in df.columns), None)

    if selected_column:
        smiles = df[[selected_column]].dropna()
    else:
        raise ValueError("SMILES, Smiles, smiles 컬럼이 존재하지 않습니다.")
    
    # 🔹 0~2047 컬럼 중 하나라도 NaN이면 삭제
    fingerprint_columns = [f"fp_{i+1}" for i in range(1, 2048)]  # 0~2047까지 컬럼 리스트
    _cal_result = process_smiles_file_ecfp6(smiles, selected_column)
    cal_result = _cal_result.dropna(subset=fingerprint_columns)
    before_drop, after_drop = len(_cal_result), len(cal_result)
    print(cal_result.head())
    print(f"\n🔹 총 {before_drop}개 중 {after_drop}개가 남음. ({before_drop - after_drop} 삭제 됨)")
    
    # 좌표 만들기
    # umap_df = umap_coordinate(cal_result, selected_column, cnt)
    umap_df = umap_coordinate(_cal_result, selected_column, cnt)
    smiles_umap_df = pd.concat([df, umap_df], axis=1)

    tsne_df = tsne_coordinate(cal_result, selected_column, cnt)
    smiles_tsne_df = pd.concat([df, tsne_df], axis=1)

    cal_df = pd.concat([df.drop(columns=selected_column), cal_result], axis=1)
    print(cal_df.head())
    cal_df.to_csv(f"{out_csv_file_path}/Rdkit_calculate_result.csv", index=False)
    
    smiles_umap_df.to_csv(f"{out_csv_file_path}/umap_result.csv", index=False)
    smiles_tsne_df.to_csv(f"{out_csv_file_path}/tSNE_result.csv", index=False)
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(
        x=smiles_umap_df['UMAP1'], 
        y=smiles_umap_df['UMAP2'], 
        hue=smiles_umap_df['tool'], 
        palette='tab10', 
        alpha=0.7, 
        s=50
    )
    plt.xlabel("UMAP 1")
    plt.ylabel("UMAP 2")
    plt.title("UMAP Scatter Plot (Colored by tool)")
    plt.legend(title="tool", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)

    # ✅ 이미지 파일로 저장
    plt.savefig(f"{out_csv_file_path}/umap.png", dpi=300, bbox_inches='tight')
    
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=smiles_tsne_df['tSNE1'], y=smiles_tsne_df['tSNE2'], hue=smiles_tsne_df['tool'], palette='tab10', alpha=0.7, s=50)
    plt.xlabel("tSNE 1")
    plt.ylabel("tSNE 2")
    plt.title("tSNE Scatter Plot (Colored by tool)")
    plt.legend(title="tool", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True)
    
    # ✅ 이미지 파일로 저장
    plt.savefig(f"{out_csv_file_path}/tsne.png", dpi=300, bbox_inches='tight')