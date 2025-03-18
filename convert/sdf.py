import sys
import logging
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import SDWriter, AllChem

logging.basicConfig(
    level=logging.INFO,
    format="%(levelname)s: %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)]  # Ensures stdout is used for INFO logs
)

def sdf(smiles: str, directory_path: str = ".", overwrite: bool = False) -> str:
    """Generates an SDF file from a SMILES string."""
    destination = Path(directory_path) / f"{smiles}.sdf"

    if not overwrite and destination.exists():
        logging.info(f"SDF already exists: {destination}")
        return str(destination)

    try:
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        mol.SetProp("SMILES", smiles)

        breakpoint()
        with SDWriter(str(destination)) as writer:
            writer.write(mol)

        logging.info(f"SDF file saved: {destination}")
        return str(destination)

    except Exception as e:
        logging.error(f"Error generating SDF for {smiles}: {e}")
        return ""

if __name__ == "__main__":
    if len(sys.argv) < 2:
        logging.error("Missing SMILES argument.")
        sys.exit(1)

    smiles_arg = sys.argv[1]
    directory_arg = sys.argv[2] if len(sys.argv) > 2 else "."
    overwrite_arg = sys.argv[3].lower() == "true" if len(sys.argv) > 3 else False
    sdf(smiles_arg, directory_arg, overwrite_arg)
