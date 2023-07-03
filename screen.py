#!/usr/bin/env python3

import typing
import subprocess
import shutil
import argparse
import logging
from pathlib import Path
import time
from functools import partial

from Bio.PDB import PDBList

import filing
import processing
import drive
import log
import chem
import rcsb_pdb


receptors_of_pdbs = {"3rzy": "fabp4", "4lkp": "fabp5"}
ref_complexes = {"fabp4": "2nnq", "fabp5": "5hz5"}
ref_ligands = {"2nnq": "T4B", "5hz5": "65X"}


def main():
    log.setup(__file__)
    args = parsed_args()
    if args.write_compounds:
        chem.compounds()
    compounds = args.compounds
    if args.default_compounds:
        compounds |= {compound_path.stem for compound_path in chem.compound_dir().iterdir()}
    compounds -= args.compound_excludes
    receptor_dir = local_data_dir() / "receptors"
    filing.mkdirs(receptor_dir, clean=True)
    for pdb_id in args.receptors_pdb_ids:
        rcsb_pdb.write_pdb(pdb_id, receptor_dir)
    rcsb_pdb.apoproteins(args.receptors, receptor_dir)
    if args.redos is not None:
        compounds |= {Path(path).stem.split('_')[0] for path in args.redos}
    if args.google_drive_destination is not None: 
        observer = processing.call_upon_file_addition(Path(args.data_dir_path) / args.dock_dir, partial(drive.Drive.write, drive_folder=Path(args.google_drive_destination)))
    docker = Docker(
        compounds=compounds, receptor_dir=receptor_dir,
        gnina_dir=args.dock_dir, default_compounds=args.default_compounds,
        sec_limit=args.sec_limit, gnina_options=args.gnina_options)
    try:
        docker.dock_batch()
    finally:
        if args.google_drive_destination is not None:
            observer.stop()
            observer.join()


def parsed_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-c", "--compounds", default=set(), nargs="*", type=str)
    parser.add_argument("-f", "--default-compounds",  action='store_true')
    parser.add_argument("-w", "--write-compounds", help="Overwrites existing compounds in the data_dir", action="store_true")
    parser.add_argument("-e", "--compound-excludes", help="Excludes from all", default={"Insulin"}, nargs="*", type=str)
    parser.add_argument("-i", "--receptors-pdb-ids", default = {"3rzy", "4lkp"}, nargs="*", type=str)
    parser.add_argument("-r", "--receptors", default=set(), help="Searches pdb database", nargs="*", type=str)
    parser.add_argument("-o", "--dock-dir", default="docks")
    parser.add_argument("-d", "--data-dir-path", default=local_data_dir())
    parser.add_argument("-s", "--sec-limit", default=float('inf'), type=float)
    parser.add_argument("-u", "--google-drive-destination", type=str)
    parser.add_argument("-a", "--redos", nargs="*", type=str)
    parser.add_argument("-g", "--gnina-options", nargs=argparse.REMAINDER, type=str, help="Override or include additional gnina option")
    return parser.parse_args()


def local_data_dir():
    return filing.local_data_dir(__file__)

class Docker:
    def __init__(
            self, compounds, data_dir=local_data_dir(),
            receptor_dir="receptors", compound_dir="compounds",
            gnina_dir="docks", default_compounds=False,
            sec_limit=float("inf"), gnina_options=[])
        self.compounds = compounds
        self.data_dir = data_dir
        self.receptor_dir = receptor_dir
        self.compound_dir = compound_dir
        self.gnina_dir = gnina_dir
        self.data_paths = {dir_name: data_dir / dir_name for dir_name in [receptor_dir, compound_dir, gnina_dir]}
        self.default_compounds = default_compounds
        self.sec_limit = sec_limit
        self.gnina_options = gnina_options

    def dock_batch(self):
        for receptor ink sorted(data_paths[receptor_dir].iterdir()):
            receptors_dock_dir: path = filing.mkdirs(data_paths[gnina_dir] / receptor.stem)
            for compound in data_paths[compound_dir].iterdir():
                if filing.matches_ignoring_spaces_and_line_symbols(compound.stem, compounds):
                    dock(receptor, compound, receptors_dock_dir, sec_limit, gnina_options)
        return data_paths[gnina_dir]


    def dock(self, receptor: path, ligand: path, destination_dir: path=path.cwd()):
        dock_result = {ext: (destination_dir / (ligand.stem + "_" + receptor.stem)).with_suffix("." + ext) for ext in ("sdf", "txt")}
        if not (dock_result["sdf"].exists() and dock_result["txt"].exists()):
            try:
                with open(dock_result["txt"], "a") as dock_txt_file:
                    logging.info(f"Starting {receptor} and {ligand}")
                    docking = subprocess.Popen(gnina_command(receptor, ligand, dock_result["sdf"], gnina_options), stdout=dock_txt_file, stderr=dock_txt_file)
                    processing.limit(docking, sec_limit, dock_txt_file)
                chem.transfer_smiles_attribute(dock_result["sdf"])
            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing {receptor} and {ligand}: {e}")
        return dock


    def gnina_command(self, receptor, ligand, gnina_result):
        gnina_path = shutil.which("gnina")
        if not gnina_path:
            gnina_path = "/home/m/cure/gnina"
        potential_autobox_ligand, receptor = autobox_ligand(receptor.stem)
        return [
            gnina_path,
            "-r", receptor,
            "-l", ligand,
            "--autobox_ligand", potential_autobox_ligand if potential_autobox_ligand is not None else receptor,
            "--autobox_extend", "1",
            "-o", gnina_result,
            "--seed", "0",
            "--exhaustiveness", "64",
        ] + self.gnina_options


    def autobox_ligand(receptor_name: str) -> dict:
        ref_complex = ref_complexes[receptors_of_pdbs[receptor_name]]
        ref_ligand = ref_ligands[ref_complex]
        complexes: Path = local_data_dir() / "complexes"
        ligand_path = (complexes / f"ligands/{ref_ligand}").with_suffix(".pdb")
        if not ligand_path.exists():
            chem.extract_ligand_from_chain(rcsb_pdb.write_pdb(ref_complex, path=complexes), "A", ref_ligand, ligand_path)
        return ligand_path


if __name__ == "__main__":
    main()

