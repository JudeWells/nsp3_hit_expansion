import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
import matplotlib.pyplot as plt

"""
Created by Jude Wells 2024-01-04
loads all molecules from Enamine Real that were similarity matched
using Infinisee "analog hunter" or "scaffold hopper" methods
search was carried out on all 10 hit compounds
duplicate molecules are left in the dataframe.


'smiles:\nid: 
name: 
result-rank: Infinisee similarity rank (1 to 500)
similarity: Infinisee similarity score to parent molecule
query-name: CACHE verified hit mol to generate similars
query-smiles: 
space: 
reaction-name: 
reagent1-name: 
reagent1-smiles: 
reagent2-name: 
reagent2-smiles: 
logp: 
mw: 
tpsa: 
parent_mol: CACHE verified hit mol to generate similars
method: {hits_analog_hunter, hits_scaffold_hopper}
parent_log10_kd_uM': KD (MicroMol) of the parent molecule used to generate the similars

"""

def display_mol(mol, title=''):
    img = Draw.MolToImage(mol)
    plt.imshow(img)
    plt.axis('off')
    plt.title(title)
    plt.show()


def contains_carboxylic_acid_or_carboxylate(smiles):
    mol = Chem.MolFromSmiles(smiles)
    carboxylate_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')
    # display_mol(carboxylate_pattern, title='carboxylate_pattern')
    return mol.HasSubstructMatch(carboxylate_pattern)

def filter_by_substructure(df):
    patt = Chem.MolFromSmarts("[c]1[n][c][n][c]2[nH][c]4[a][a][a][a][a]4[a]12")
    display_mol(patt, title='scaffold_pattern')
    drop_indices = []
    for i, row in df.iterrows():
        mol = row['mol']
        if not mol.HasSubstructMatch(patt):
            # display_mol(mol)
            drop_indices.append(i)
    print(f"dropping {len(drop_indices)} molecules due to lack of scaffold")
    df.drop(drop_indices, inplace=True)
    return df

def remove_carboxylates(df):
    carboxylate_pattern = Chem.MolFromSmarts('[CX3](=O)[OX1H0-,OX2H1]')
    # display_mol(carboxylate_pattern)
    drop_indices = []
    for i, row in df.iterrows():
        mol = row['mol']
        if mol.HasSubstructMatch(carboxylate_pattern):
            # display_mol(mol)
            drop_indices.append(i)
    print(f"dropping {len(drop_indices)} molecules due to carboxylate")
    df.drop(drop_indices, inplace=True)
    return df

if __name__=="__main__":
    df = pd.read_csv('cache3_analog_hunters_and_scaffold_hoppers_from_10_hits.csv')
    print(f"Total Number of molecules (all matches including duplicates) {len(df)}")
    df['mol'] = [Chem.MolFromSmiles(smile) for smile in df['smiles']]
    df = filter_by_substructure(df)
    df = remove_carboxylates(df)
    print(f"Total Number of molecules (all matches) after carboxyl+substructure filtering {len(df)}")
    print(f"Total Number of unique molecules after carboxyl+substructure filtering {len(df['smiles'].unique())}")