from pathlib import 
import logging
from typing import Dict, List

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import SDWriter
from Bio.PDB import PDBIO, PDBParser, Select, Structure, Model, Chain, Atom
from Bio.PDB.Residue import Residue

import spreadsheet
import filing


def local_data_dir():
    return filing.local_data_dir(__file__)


def compounds_dataframe() -> pd.DataFrame:
    smiless_with_names_range = "A:B"
    names_smiless_value_range = spreadsheet.get_values(
        "1VRl8dghA5t1JiLEa47OWW-hr-Qbvn9yksB4Hrga19Ck", smiless_with_names_range
    )
    compounds_spreadsheet = names_smiless_value_range["values"]
    compounds_dataframe = pd.DataFrame(compounds_spreadsheet[4:], columns=compounds_spreadsheet[3])
    compounds_dataframe["Name"] = compounds_dataframe["Name"].replace('', np.nan) 
    compounds_dataframe["Name"] = compounds_dataframe["Name"].fillna("index_" + compounds_dataframe.index.to_series().astype(str))
    return compounds_dataframe


def compound_dir():
    return (__file__).parent / "data/compounds"


def compounds(directory_path=compound_dir(), excludes=[], overwrite=False):
    for _, name, smiles in compounds_dataframe().itertuples():
        if name not in excludes:
            sdf(name.replace(' ', '_'), smiles, overwrite, directory_path)
    return directory_path


def sdf(name, smiles, overwrite=False, directory_path=compound_dir()):
    destination = directory_path / (name).with_suffix(".sdf")
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


class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.get_full_id()[3][0] == ' '


class LigandSelect(Select):
    def __init__(self, chain_id, residue_name):
        self.chain_id = chain_id
        self.residue_name = residue_name

    def accept_chain(self, chain):
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue):
        return residue.get_resname() == self.residue_name


def extract_from_structure(pdb_path: Path, output_file: Path, select: Select):
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    return write_structure(structure, output_file, select)


def extract_residue(pdb_path: Path, identifier: str, output_file: Path) -> Dict[str, List[Residue]]:
    parser = PDBParser()
    structure = parser.get_structure(pdb_path.stem, str(pdb_path))
    for residue in structure.get_residues():
        if residue.get_resname() == identifier.upper():
            write_residue(residue, output_file)
            return residue


def write_residue(residue: Residue, output_file: ) -> Path:
    write_structure(add_recursively(*ligand_scaffold(), residue), output_file)
    return output_file


def write_structure(structure: Structure, output_file: , select: Select=None) -> Path:
    pdb_io = PDBIO()
    pdb_io.set_structure(structure)
    pdb_io.save(str(output_file), select)
    return output_file


def ligand_scaffold():
    ligand_structure = Structure.Structure("Ligand Structure")
    model = Model.Model(0)
    chain = Chain.Chain("A")
    return (ligand_structure, model, chain)


def add_recursively(structure, model, chain, residue, atom):
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


def pdb_partition(pdb_path: , identifier: str):
    pdb_partition_path = (pdb_path.parent / identifier).with_suffix(".pdb")
    with open(pdb_partition_path, "w") as pdb_path_file:
        subprocess.run(["grep", identifier, pdb_path], stdout=pdb_path_file)
    return pdb_partition_path


def transfer_smiles_attributes(): 
    docks_dir = local_data_dir() / "docks_drug_autobox" / "fabp4_apo_3RZY"
    for dock_path in docks_dir.iterdir():
        try:
            transfer_smiles_attribute(dock_path)
        except:
            continue


def transfer_smiles_attribute(dock: ):
    source_supplier = Chem.SDMolSupplier(str(compound_of_dock(dock)))
    for mol in source_supplier:
        if mol is not None:
            smiles = mol.GetProp('SMILES')
    dock_suppliers_mols = [mol for mol in Chem.SDMolSupplier(str(dock)) if mol is not None]
    with Chem.SDWriter(str(dock)) as writer:
        for mol in dock_suppliers_mols:
            mol.SetProp('SMILES', smiles)
            writer.write(mol)


def compound_of_dock(dock: ) -> Path:
    compound = '_'.join(dock.stem.split('_')[:-1])
    return (local_data_dir() / "compounds" / compound).with_suffix(".sdf")

