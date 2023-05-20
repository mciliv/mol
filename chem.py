from pathlib import Path
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

import log

log.setup(__file__)

from pathlib import Path

import util
import chem
import spreadsheetor


def compound_names_smiless():
    spreadsheet = spreadsheetor.get_spreadsheet()
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheetor.get_values(spreadsheetor.SPREADSHEET_ID, smiless_with_names_range)
    names_smiless = names_smiless_value_range['values']
    names_smiless = list(zip(*names_smiless))
    data_start = 4
    return names_smiless[data_start:]


def compound_dir():
    return Path(__file__).parent / "data/compounds"

def compounds(directory_path=compound_dir()):
    for compound in compound_names_smiless():
        chem.sdf(*compound, directory_path=directory_path)
    return directory_path


def sdf(name, smiles, overwrite=False, directory_path=compound_dir()):
    destination = directory_path / Path(name).with_suffix(".sdf")
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




