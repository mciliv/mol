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
import log
import chem
import rcsb_pdb


def main():
    log.setup(__file__)
    args = parsed_args()
    print(args.dock_dir)
    print(args.dock_dir)
    exit()
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
    dock_dir = args.dock_dir
    if args.redos is not None:
        compounds |= {Path(path).stem.split('_')[0] for path in args.redos}
    #observer = processing.start_watchdog(dock_dir, partial(drive.Drive.write(
    try:
        dock_batch(
            compounds=compounds, receptor_dir=receptor_dir, gnina_dir=args.dock_dir, default_compounds=args.default_compounds, sec_limit=args.sec_limit
        )
    finally:
        observer.stop()
        observer.join()


def parsed_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--compounds", default=set(), nargs="*", type=str)
    parser.add_argument("-w", "--write-compounds", help="Overwrites existing compounds", action="store_true")
    parser.add_argument("--default-compounds",  action='store_true')
    parser.add_argument("-e", "--compound-excludes", help="Excludes from all", default={"Insulin"}, nargs="*", type=str)
    parser.add_argument("-i", "--receptors-pdb-ids", default = {"3rzy", "4lkp"}, nargs="*", type=str)
    parser.add_argument("-r", "--receptors", default=set(), help="Searches pdb database", nargs="*", type=str)
    parser.add_argument("-o", "--dock-dir", default="docks")
    parser.add_argument("-d", "--data-dir-path", default=local_data_dir())
    parser.add_argument("-s", "--sec-limit", default=float('inf'), type=float)
    parser.add_argument("-u", "--google-drive-destination", type=str)
    parser.add_argument("-a", "--redos", nargs=argparse.REMAINDER, type=str)
    return parser.parse_args()


def local_data_dir():
    return filing.local_data_dir(__file__)


def dock_batch(compounds, receptor_dir="receptors", compound_dir="compounds", gnina_dir="docks", default_compounds=False, sec_limit=float("inf")):
    data_paths = {dir: local_data_dir() / dir for dir in [receptor_dir, compound_dir, gnina_dir]}
    for receptor in sorted(data_paths[receptor_dir].iterdir()):
        receptors_dock_dir: Path = filing.mkdirs(data_paths[gnina_dir] / receptor.stem)
        for compound in data_paths[compound_dir].iterdir():
            if filing.matches_ignoring_spaces_and_line_symbols(compound.stem, compounds):
                dock(receptor, compound, receptors_dock_dir, sec_limit)
    return data_paths[gnina_dir]


def dock(receptor: Path, ligand: Path, destination_dir: Path=Path.cwd(), sec_limit=float('inf')):
    dock_result = {ext: (destination_dir / (ligand.stem + "_" + receptor.stem)).with_suffix("." + ext) for ext in ("sdf", "txt")}
    if not (dock_result["sdf"].exists() and dock_result["txt"].exists()):
        try:
            with open(dock_result["txt"], "a") as dock_txt_file:
                logging.info(f"Starting {receptor} and {ligand}")
                breakpoint()
                docking = subprocess.Popen(
                        gnina_command(receptor, ligand, dock_result["sdf"]),
                        stdout=dock_txt_file,
                        stderr=dock_txt_file,
                        )
                processing.limit(docking, sec_limit, dock_txt_file)
            chem.transfer_smiles_attribute(dock_result["sdf"])
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {receptor} and {ligand}: {e}")
    return dock


def gnina_command(receptor, ligand, gnina_result):
    gnina_path = shutil.which("gnina")
    if not gnina_path:
        gnina_path = "/home/m/cure/gnina"
    potential_autobox_ligand = autobox_ligand(receptor.stem)
    return [
        gnina_path,
        "-r",
        receptor,
        "-l",
        ligand,
        "--autobox_ligand",
        potential_autobox_ligand if potential_autobox_ligand is not None else receptor,
        "--autobox_extend",
        "1",
        "-o",
        gnina_result,
        "--seed",
        "0",
        "--exhaustiveness",
        "64",
    ]


def autobox_ligand(receptor_name: str) -> dict:
    receptors_of_pdbs = {"3rzy": "fabp4", "4lkp": "fabp5"}
    ref_complexes = {"fabp4": "2nnq", "fabp5": "5hz5"}
    ref_ligands = {"2nnq": "T4B", "5hz5": "A65X"}
    ref_complex = ref_complexes[receptors_of_pdbs[receptor_name]]
    ref_ligand = ref_ligands[ref_complex]
    complexes: Path = local_data_dir() / "complexes"
    ligand_path = (complexes / f"ligands/{ref_ligand}").with_suffix(".pdb")
    if not ligand_path.exists():
        chem.extract_ligand(rcsb_pdb.write_pdb(ref_complex, path=complexes), ref_ligand, ligand_path)
    return ligand_path


if __name__ == "__main__":
    main()

