import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolTransforms
import pandas as pd


# Function to calculate dihedral angle
def calculate_dihedral(mol, atom_indices):
    conf = mol.GetConformer()
    dihedral_angle = rdMolTransforms.GetDihedralDeg(conf, *atom_indices)
    return dihedral_angle


# Sample molecules
smiles_list = ["CCO", "CCCC", "CCN", "CCOCC"]
molecules = [Chem.MolFromSmiles(smiles) for smiles in smiles_list]

# Generate 3D conformers
for mol in molecules:
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)

# Collect data
data = []
for mol in molecules:
    atom_coords = mol.GetConformer().GetPositions()
    for i in range(mol.GetNumAtoms()):
        atom_symbol = mol.GetAtomWithIdx(i).GetSymbol()
        coords = atom_coords[i]
        for j in range(i + 1, mol.GetNumAtoms()):
            for k in range(j + 1, mol.GetNumAtoms()):
                for l in range(k + 1, mol.GetNumAtoms()):
                    try:
                        dihedral = calculate_dihedral(mol, [i, j, k, l])
                        data.append(
                            {
                                "SMILES": Chem.MolToSmiles(mol),
                                "Atom 1": i,
                                "Atom 2": j,
                                "Atom 3": k,
                                "Atom 4": l,
                                "Dihedral Angle (deg)": dihedral,
                                "Atom 1 Coord": atom_coords[i],
                                "Atom 2 Coord": atom_coords[j],
                                "Atom 3 Coord": atom_coords[k],
                                "Atom 4 Coord": atom_coords[l],
                            }
                        )
                    except:
                        continue

# Convert to DataFrame
df = pd.DataFrame(data)

# Show the data
print(df.head())

# Optionally, save to a CSV file
df.to_csv("./data/molecule_dihedrals.csv", index=False)
