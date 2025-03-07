import sys
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem

logging.basicConfig(level=logging.INFO)

def sdf(smiles: str, directory_path: str = ".", overwrite: bool = False):
    """Generates an SDF file from a SMILES string."""
    destination = Path(directory_path) / f"{smiles}.sdf"

    if not overwrite and destination.exists():
        return str(destination)

    try:
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        Chem.RemoveHs(mol)
        mol.SetProp("SMILES", smiles)

        with SDWriter(str(destination)) as writer:
            writer.write(mol)

        logging.info(f"SDF file saved: {destination}")

    except Exception as e:
        logging.error(f"Error generating SDF for {smiles}: {e}")

    return destination


def is_valid_sdf_with_molecule(filepath):
    try:
        supplier = Chem.SDMolSupplier(str(filepath))
        for mol in supplier:
            if mol:
                return True
    finally:
        return False


if __name__ == "__main__":
    smiles_arg = sys.argv[1]
    directory_arg = sys.argv[2] if len(sys.argv) > 2 else "."
    overwrite_arg = sys.argv[3].lower() == "true" if len(sys.argv) > 3 else False
    sdf(smiles_arg, directory_arg, overwrite_arg)

