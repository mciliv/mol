#!/usr/bin/env python3

import subprocess
import argparse
import logging
from pathlib import Path
import write
import time
from datetime import datetime

import log
from util.util import clean_directory

log.setup(__file__)


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--write-compounds', action='store_true')
    parser.add_argument('-n', '--no-dock', action='store_true')
    parser.add_argument('-o', '--output-dir', default=Path(__file__).resolve().parent)
    return parser.parse_args()


def run_gnina_on_permutations(ligand_dir, receptor_pdbs, output_dir):
    if not ligand_dir.exists():
        logging.error("Candidates directory not found.")
        return

    for receptor, receptor_pdb in receptor_pdbs.items():
        if not receptor_pdb.exists():
            logging.error(receptor_pdb + " not found.")
            return
        receptors_dock_dir = output_dir / receptor
        receptors_dock_dir.mkdir(parents=True, exist_ok=True)

        for ligand_sdf in ligand_dir.iterdir():
            dock_file_no_ext = receptors_dock_dir / (ligand_sdf.stem + "_" + receptor_pdb.stem)
            dock_sdf = dock_file_no_ext.with_suffix(".sdf")
            dock_txt = dock_file_no_ext.with_suffix(".txt")

            command = ["./gnina", "-r", receptor_pdb,
                       "-l", ligand_sdf,
                       "--autobox_ligand", ligand_sdf,
                       "-o", dock_sdf,
                       "--seed", "0"]

            if not (dock_sdf.exists() and dock_txt.exists()):
                try:
                    with open(dock_txt, 'a') as dock_txt_file:
                        logging.info(f"Starting {receptor} and {ligand_sdf}")
                        docking = subprocess.Popen(command, stdout=dock_txt_file, stderr=dock_txt_file)
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
                    logging.error(f"Error processing {receptor} and {ligand}: {e}")


if __name__ == "__main__":
    args = args()
    output_dir = Path(args.output_dir)
    ligand_dir = output_dir / "ligands"
    gene_pdbs = write.gene_pdbs(['fabp4', 'fabp5'], output_dir / "receptors")
    if args.write_compounds:
        write.compounds(ligand_dir)
    if not args.no_dock:
        run_gnina_on_permutations(ligand_dir, gene_pdbs, output_dir / f"docks_{datetime.now().strftime('%Y-%m-%d_%H-%M')}")
