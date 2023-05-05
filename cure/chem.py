from pathlib import Path
import logging

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter

import log

log.setup(__file__)


def sdf(name, smiles, path='./'):
    try:
        mol = Chem.MolFromSmiles(smiles)
        Chem.Kekulize(mol)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        Chem.RemoveHs(mol)
        with SDWriter(path / Path(name).with_suffix(".sdf")) as writer:
            writer.write(mol)
    except Exception as e:
        logging.error("For " + name + ":\n" + str(e))

