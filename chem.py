from pathlib import Path
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

import log

log.setup(__file__)

from pathlib import Path

from util.util import clean_directory
import chem, spreadsheetor


def compound_names_smiles():
    spreadsheet = spreadsheetor.get_spreadsheet()
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheetor.get_values(spreadsheetor.SPREADSHEET_ID, smiless_with_names_range)
    names_smiless = names_smiless_value_range['values']
    names_smiless = list(zip(*names_smiless))
    data_start = 4
    return names_smiless[data_start:]


def compounds(path=Path(__file__).parent / "compounds"):
    for compound in compound_names_smiles():
        chem.sdf(*compound, path=path)
    return path


def sdf(name, smiles, overwrite=False, path=Path(__file__).parent):
    destination = path / Path(name).with_suffix(".sdf")
    if overwrite or not destination.exists():
        try:
            mol = Chem.MolFromSmiles(smiles)
            Chem.Kekulize(mol)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            Chem.RemoveHs(mol)
            with open(destination, 'w') as file:
                with SDWriter(file) as writer:
                    writer.write(mol)
        except Exception as e:
            logging.error(f"\n{name}\n{e}\n")
    return destination




