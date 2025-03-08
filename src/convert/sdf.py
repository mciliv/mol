import sys
import logging
import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem

logging.basicConfig(level=logging.INFO)


def sdf(smiles: str, directory_path: str = ".", overwrite: bool = False) -> str:
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

    return str(destination)


def is_valid_sdf_with_molecule(filepath):
    try:
        supplier = Chem.SDMolSupplier(str(filepath))
        for mol in supplier:
            if mol:
                return True
    finally:
        return False


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate an SDF file from a SMILES string.")
    parser.add_argument("smiles", type=str, help="SMILES string for the molecule")
    parser.add_argument("--dir", type=str, default=".", help="Directory to save the SDF file")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing file if set")

    args = parser.parse_args()

    sdf(args.smiles, args.dir, args.overwrite)
