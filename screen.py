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

log.setup(__file__)


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--redock', type=str)
    parser.add_argument('-c', '--write-compounds', action='store_true')
    parser.add_argument('-a', '--dock-all', action='store_true', help='Docks all combos')
    parser.add_argument('-o', '--output-dir', default=Path(__file__).resolve().parent / "data")
    return parser.parse_args()


def data_dir():
    return util.data_dir(__file__)

def gnina_command(receptor_pdb, ligand_sdf, output_sdf):
    gnina_path = shutil.which("gnina")
    if not gnina_path:
        gnina_path = "/home/m/cure/gnina"
    return [gnina_path, "-r", receptor_pdb,
               "-l", ligand_sdf,
               "--autobox_ligand", ligand_sdf,
               "-o", output_sdf,
               "--seed", "0"]


def dock(receptor: Path, ligand: Path, output_stem: Path, sec_limit=60):
    dock = {file_type: output_stem.with_suffix("." + file_type) for file_type in ("sdf", "txt")}
    if not (dock["sdf"].exists() and dock["txt"].exists()):
        try:
            with open(dock["txt"], 'a') as dock_txt_file:
                logging.info(f"Starting {receptor} and {ligand}")
                docking = subprocess.Popen(gnina_command(receptor, ligand, dock["sdf"]), stdout=dock_txt_file,
                                            stderr=dock_txt_file)
                stopped = False
                for duration in range(1, sec_limit + 1):
                    time.sleep(1)
                    poll = docking.poll()
                    if poll is not None:
                        dock_txt_file.write(f"Ran for {duration} sec; process poll value is {poll}")
                        stopped = True
                        logging.info(f"Finished: {output_stem}")
                        break
                if not stopped:
                    dock_txt_file.write(f"Terminated with {sec_limit} sec limit")
        except subprocess.CalledProcessError as e:
            logging.error(f"Error processing {receptor} and {ligand}: {e}")
    return dock


def dock_all(receptor_dir='receptors', compound_dir='compounds', output_dir='docks'):
    data_paths = {dir: data_dir() / dir for dir in [receptor_dir, compound_dir, output_dir]}
    for receptor in sorted(data_paths[receptor_dir].iterdir()):
        receptors_dock_dir = util.mkdirs(data_paths[output_dir] / receptor.stem)

        for ligand_sdf in data_paths[compound_dir].iterdir():
            output_stem_path = receptors_dock_dir / (ligand_sdf.stem + "_" + receptor.stem)
            dock(receptor, ligand_sdf, output_stem_path)


def pdb_partition(pdb_path: Path, identifier: str):
    pdb_partition_path = (pdb_path.parent / identifier).with_suffix('.pdb')
    with open(pdb_partition_path, "w") as pdb_path_file:
        subprocess.run(["grep", identifier, pdb_path], stdout=pdb_path_file)
    return pdb_partition_path


if __name__ == "__main__":
    args = args()
    if args.write_compounds:
        chem.compounds()
    if args.dock_all:
        dock_all(rcsb_pdb.apoprotein(['fabp4', 'fabp5']).stem)



