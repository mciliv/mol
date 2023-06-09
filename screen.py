#!/usr/bin/env python3

import typing
import subprocess
import shutil
import argparse
import logging
from pathlib import Path
import time
from datetime import datetime

import filing
import log
import chem
import rcsb_pdb


def parsed_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--compounds", type=str, nargs="*", default=[])
    parser.add_argument("-w", "--write-compounds", action="store_true")
    parser.add_argument("--all", action='store_true')
    parser.add_argument("-e", "--exclude-compounds", type=str, nargs="*", default=[])
    parser.add_argument("-r", "--receptors", type=str, nargs="*", default=["fabp4", "fabp5"])
    parser.add_argument(
        "-o",
        "--dock-dir",
        default=f"docks_ligand_autoboxed_{datetime.now().isoformat()}",
    )
    parser.add_argument(
        "-d", "--data-dir", default=local_data_dir()
    )
    parser.add_argument(
        "-a", "--redos", type=str, nargs=argparse.REMAINDER
    )
    return parser.parse_args()


def local_data_dir():
    return filing.local_data_dir(__file__)


def dock_batch(compounds, receptor_dir="receptors", compound_dir="compounds", gnina_dir="docks", dock_all=False):
    data_paths = {
        dir: local_data_dir() / dir for dir in [receptor_dir, compound_dir, gnina_dir]
        }
    for receptor in sorted(data_paths[receptor_dir].iterdir()):
        receptors_dock_dir: Path = filing.mkdirs(data_paths[gnina_dir] / receptor.stem)

        for compound in data_paths[compound_dir].iterdir():
            if not dock_all and not filing.matches_ignoring_spaces_and_line_symbols(compound.stem, compounds):
                continue
            dock(receptor, compound, receptors_dock_dir)


def dock(receptor: Path, ligand: Path, destination_dir: Path=Path.cwd(), sec_limit=float('inf')):
    dock_result = \
            {file_type: receptors_dock_dir / (
                    ligand.stem + "_" + receptor.stem
                    ).with_suffix("." + file_type)
    
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
                sec = 0
                while sec < sec_limit:
                    time.sleep(1)
                    poll = docking.poll()
                    if poll is not None:
                        dock_txt_file.write(
                                f"Ran for {sec} sec; process poll value is {poll}"
                                )
                        stopped = True
                        logging.info(f"Finished: {gnina_result_stem}")
                        break
                    sec += 1
                if not stopped:
                    dock_txt_file.write(f"Terminated with {sec_limit} sec limit")
            chem.transfer_smiles_attribute(ligand, dock_result["sdf"])
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
        "-o",
        gnina_result,
        "--seed",
        "0",
        "--exhaustiveness",
        "64",
    ]


def autobox_ligand(receptor_name: str) -> dict:
    receptors_of_pdbs = {"3RZY": "fabp4", "4LKP": "fabp5"}
    ref_complexes = {"fabp4": "2NNQ", "fabp5": "5HZ5"}
    ref_ligands = {"2NNQ": "T4B", "5HZ5": "65X"}
    ref_complex = ref_complexes[receptors_of_pdbs[receptor_name]]
    ref_ligand = ref_ligands[ref_complex]
    complexes: Path = local_data_dir() / "complexes"
    ligand_path = (complexes / f"ligands/{ref_ligand}").with_suffix(".pdb")
    if not ligand_path.exists():
        chem.extract_ligand(rcsb_pdb.write_pdb(ref_complex, path=complexes), ref_ligand, ligand_path)
    return ligand_path


if __name__ == "__main__":
    log.setup(__file__)
    args = parsed_args()
    if args.write_compounds:
        chem.compounds()
    compounds = args.compounds
    if args.redos is not None:
        compounds += [Path(path).stem.split('_')[0] for path in args.redos]
    dock_batch(
        compounds=compounds, receptor_dir=rcsb_pdb.apoproteins(args.receptors), gnina_dir=args.dock_dir, dock_all=args.dock_all
    )


