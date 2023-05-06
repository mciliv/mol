#!/usr/bin/env python3

import subprocess
import argparse
import logging
from pathlib import Path
import write
import time
from datetime import datetime

import log
import rcsb_pdb

log.setup(__file__)


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--redock', type=str)
    parser.add_argument('-c', '--write-compounds', action='store_true')
    parser.add_argument('-a', '--dock-all', action='store_true', description='Docks all in receptor and ligands folder')
    parser.add_argument('-o', '--output-dir', default=Path(__file__).resolve().parent)
    return parser.parse_args()


def gnina_command(dock_sdf, ligand_sdf, receptor_pdb):
    return ["./gnina", "-r", receptor_pdb,
               "-l", ligand_sdf,
               "--autobox_ligand", ligand_sdf,
               "-o", dock_sdf,
               "--seed", "0"]


def dock_all(receptor_dir=Path('./receptors/'), compound_dir=Path('./compounds/'), output_dir=Path(f"./docks_{datetime.now().strftime('%Y-%m-%d_%H-%M')}/")):
    if not compound_dir.exists():
        logging.error("Candidates directory not found.")
        return

    for receptor in receptor_dir.iterdir():
        receptors_dock_dir = output_dir / receptor
        receptors_dock_dir.mkdir(parents=True, exist_ok=True)

        for ligand_sdf in compound_dir.iterdir():
            dock_file_no_ext = receptors_dock_dir / (ligand_sdf.stem + "_" + receptor.stem)
            dock_sdf = dock_file_no_ext.with_suffix(".sdf")
            dock_txt = dock_file_no_ext.with_suffix(".txt")

            if not (dock_sdf.exists() and dock_txt.exists()):
                try:
                    with open(dock_txt, 'a') as dock_txt_file:
                        logging.info(f"Starting {receptor} and {ligand_sdf}")
                        docking = subprocess.Popen(gnina_command(dock_sdf, ligand_sdf, receptor), stdout=dock_txt_file, stderr=dock_txt_file)
                        limit = 60
                        stopped = False
                        for duration in range(1, limit + 1):
                            time.sleep(1)
                            poll = docking.poll()
                            if poll is not None:
                                dock_txt_file.write(f"Ran for {duration} sec; process poll value is {poll}")
                                stopped = True
                                logging.info(f"Finished: {dock_file_no_ext}")
                                break
                        if not stopped:
                            dock_txt_file.write(f"Terminated with {limit} sec limit")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Error processing {receptor} and {compound_dir}: {e}")

def seperation_command():
    command = [
        "grep", "ATOM", output_path.stem, '>', 'rec.pdb'
    ]

def redock(pdb_id):
    output_path = rcsb_pdb.write_pdb(pdb_id, prepended_label="redock")

    subprocess.run()


if __name__ == "__main__":
    args = args()
    if args.write_compounds:
        write.compounds()
    if args.dock_all:
        dock_all(rcsb_pdb.apoprotein(['fabp4', 'fabp5']))



