from pathlib import Path
import logging
from typing import Dict, List

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from Bio.PDB import PDBIO, PDBParser, Structure, Model, Chain, Atom
from Bio.PDB.Residue import Residue

import spreadsheet


def compound_names_smiless() -> pd.DataFrame:
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheet.get_values(
        "1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck", smiless_with_names_range
    )
    compounds_spreadsheet = names_smiless_value_range["values"]
    return pd.DataFrame(compounds_spreadsheet[4:], columns=compounds_spreadsheet[3])


def compound_dir():
    return Path(__file__).parent / "data/compounds"


def compounds(directory_path=compound_dir()):
    for compound in ChemSource.compound_names_smiless().iterrows():
        sdf(compound["name"].replace(' ', '_'), compound["smiles"], directory_path=directory_path)
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
            
            mol.SetProp("_Name", name) 
            mol.SetProp("SMILES", smiles) 
            
            with open(destination, "w") as file:
                with SDWriter(file) as writer:
                    writer.write(mol)
        except Exception as e:
            logging.error(f"\n{name}\n{e}\n")
    return destination


def extract_ligand(pdb_file: Path, identifier: str, output_file: Path) -> Dict[str, List[Residue]]:
    parser = PDBParser()
    structure = parser.get_structure(pdb_file.stem, str(pdb_file))
    for residue in structure.get_residues():
        if residue.get_id()[0] == f'H_{identifier.upper()}':
            write_ligand(residue, output_file)
            return residue


def write_ligand(residue: Residue, output_file: Path) -> Path:
    write_structure(add_recursively(*ligand_scaffold(), residue), output_file)
    return output_file


def write_structure(structure: Structure, output_file: Path) -> Path:
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(str(output_file))
    return output_file


def ligand_scaffold():
    ligand_structure = Structure.Structure("Ligand Structure")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    return (ligand_structure, model, chain)


def add_recursively(structure, model, chain=None, residue=None, atom=None):
    if atom is not None:
        residue.add(atom)
    chain.add(residue)
    model.add(chain)
    structure.add(model)
    return structure


def accumulate_parts(parts, super_part=Residue((" ", 0, " "), "LIG", " ")):
    for part in parts:
        super_part.add(part)
    return super_part


def pdb_partition(pdb_path: Path, identifier: str):
    pdb_partition_path = (pdb_path.parent / identifier).with_suffix(".pdb")
    with open(pdb_partition_path, "w") as pdb_path_file:
        subprocess.run(["grep", identifier, pdb_path], stdout=pdb_path_file)
    return pdb_partition_path

