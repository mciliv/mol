from rdkit import Chem
from rdkit.Chem import SDWriter

def sdf(smiles, overwrite=False, directory_path='.'):
    destination = directory_path / Path(smiles).with_suffix(".sdf")
    if overwrite or not destination.exists():
        try:
            mol = Chem.MolFromSmiles(smiles)
            Chem.Kekulize(mol)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            Chem.RemoveHs(mol)
            
            mol.SetProp("SMILES", smiles) 
            
            with open(destination, "w") as file:
                with SDWriter(file) as writer:
                    writer.write(mol)
        except Exception as e:
            logging.exception(f"\n{name}\n{e}\n")
    return destination


def is_valid_sdf_with_molecule(filepath):
    try:
        supplier = Chem.SDMolSupplier(str(filepath))
        for mol in supplier:
            if mol:
                return True
    finally:
        return False

sdf(argv[1])
