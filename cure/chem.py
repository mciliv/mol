import os
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

def sdf(name, smiles, path='./'):
    errors = os.path.join(path, "errors.txt")
    try:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        Chem.RemoveHs(mol)
        with SDWriter(os.path.join(path, name + ".sdf")) as writer:
            writer.write(mol)
    except Exception as e:
        mode = 'a' if os.path.exists(errors) else 'w'
        with open(errors, mode) as f:
            print(e)
            print(e, file=f)
