#!/usr/bin/env python3

import typing
import subprocess
import shutil
import argparse
import logging
from pathlib import Path
import time
from datetime import datetime

import util
import log
import chem
import rcsb_pdb


def parsed_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-w", "--write-compounds", action="store_true")
    parser.add_argument("-c", "--compounds", type=str, nargs="+")
    parser.add_argument("-p", "--proteins", type=str, nargs="+", default=["fabp4", "fabp5"])
    parser.add_argument(
        "-o",
        "--dock-dir",
        default=f"docks_ligand_autoboxed_{datetime.now().isoformat()}",
    )
    parser.add_argument(
        "-a", "--dock-all", action="store_true", help="Docks all combos"
    )
    parser.add_argument(
        "-d", "--data-dir", default=Path(__file__).resolve().parent / "data"
    )
    return parser.parse_args()


def local_data_dir():
    return util.local_data_dir(__file__)


def autobox_ligand(protein: str) -> dict:
    ref_complexes = {"fabp4": "2NNQ", "fabp5": "5HZ5"}
    ref_ligands = {"2NNQ": "T4B", "5HZ5": "65X"}
    complexes: Path = local_data_dir() / "complexes"
    ligand_path = (
        complexes / f"ligands/{ref_ligands[ref_complexes[protein]]}"
    ).with_suffix(".pdb")
    if not ligand_path.exists():
        ref_complex = ref_complexes[protein]
        ref_ligand = ref_ligands[ref_complex]
        chem.extract_ligand(rcsb_pdb.write_pdb(ref_complex, path=complexes), ref_ligand, ligand_path)
    return ligand_path


def gnina_command(receptor, ligand, gnina_result):
    gnina_path = shutil.which("gnina")
    if not gnina_path:
        gnina_path = "/home/m/cure/gnina"
    protein = receptor.stem[:5]
    potential_autobox_ligand = autobox_ligand(protein)
    return [
        gnina_path,
        "-r",
        receptor,
        "-l",
        ligand,
        "--autobox_ligand",
        potential_autobox_ligand if potential_autobox_ligand is not None else receptor,
        "-o",
        gnina_result,
        "--seed",
        "0",
        "--exhaustiveness",
        "64",
    ]


def dock(receptor: Path, ligand: Path, gnina_result_stem: Path, sec_limit=480):
    dock_result = {
        file_type: gnina_result_stem.with_suffix("." + file_type)
        for file_type in ("sdf", "txt")
    }
    if not (dock_result["sdf"].exists() and dock_result["txt"].exists()):
        try:
            with open(dock_result["txt"], "a") as dock_txt_file:
                logging.info(f"Starting {receptor} and {ligand}")
                docking = subprocess.Popen(
                    gnina_command(receptor, ligand, dock_result["sdf"]),
                    stdout=dock_txt_file,
                    stderr=dock_txt_file,
                )
                stopped = False
                for duration in range(1, sec_limit + 1):
                    time.sleep(1)
                    poll = docking.poll()
                    if poll is not None:
                        dock_txt_file.write(
                            f"Ran for {duration} sec; process poll value is {poll}"
                        )
                        stopped = True
                        logging.info(f"Finished: {gnina_result_stem}")
                        break
                if not stopped:
                    dock_txt_file.write(f"Terminated with {sec_limit} sec limit")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {receptor} and {ligand}: {e}")
    return dock


def dock_all(receptor_dir="receptors", compound_dir="compounds", gnina_dir="docks"):
    data_paths = {
        dir: local_data_dir() / dir for dir in [receptor_dir, compound_dir, gnina_dir]
    }
    for receptor in sorted(data_paths[receptor_dir].iterdir()):
        receptors_dock_dir: Path = util.mkdirs(data_paths[gnina_dir] / receptor.stem)

        for ligand in data_paths[compound_dir].iterdir():
            gnina_result_stem: Path = receptors_dock_dir / (
                ligand.stem + "_" + receptor.stem
            )
            dock(receptor, ligand, gnina_result_stem)


def pdb_partition(pdb_path: Path, identifier: str):
    pdb_partition_path = (pdb_path.parent / identifier).with_suffix(".pdb")
    with open(pdb_partition_path, "w") as pdb_path_file:
        subprocess.run(["grep", identifier, pdb_path], stdout=pdb_path_file)
    return pdb_partition_path


if __name__ == "__main__":
    log.setup(__file__)
    args = parsed_args()
    if args.write_compounds:
        chem.compounds()
    if args.dock_all:
        dock_all(
            receptor_dir=rcsb_pdb.apoprotein(args.proteins), gnina_dir=args.dock_dir
        )
