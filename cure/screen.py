#!/usr/bin/env python3

import subprocess
import argparse
import logging
from pathlib import Path
import write

from util.util import clean_directory

logging.basicConfig(level=logging.INFO)
file_handler = logging.FileHandler(Path.home() / 'cure/cure/logs.log')
logging.getLogger().addHandler(file_handler)


def args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir')
    return parser.parse_args()


def run_gnina_on_permutations(candidates_dir, receptor_pdbs, output_dir):
    if not candidates_dir.exists():
        logging.error("Candidates directory not found.")
        return

    clean_directory(output_dir)

    for receptor, receptor_pdb in receptor_pdbs.items():
        if not receptor_pdb.exists():
            logging.error(receptor_pdb + " not found.")
            return
        receptor_output_dir = output_dir / receptor
        clean_directory(receptor_output_dir)

        for candidate in candidates_dir.iterdir():
            candidate_sdf = candidates_dir / candidate
            output_file = candidate.stem + "_" + receptor_pdb.stem + ".sdf"
            output_path = receptor_output_dir / output_file

            command = ["../../gnina", "-r", receptor_pdb,
                       "-l", candidate_sdf,
                       "--autobox_ligand", candidate_sdf,
                       "-o", output_path,
                       "--seed 0"]
            try:
                subprocess.run(command, check=True)
                logging.info(f"Processed {receptor} and {candidate}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error processing {receptor} and {candidate}: {e}")


if __name__ == "__main__":
    project_dir = Path(args().output_dir)
    candidates_dir = project_dir / "ligands"
    receptors_dir = project_dir / "receptors"
    output_dir = project_dir / "docks"

    gene_pdbs = write.gene_pdbs(['fabp4', 'fabp5'], receptors_dir)
    run_gnina_on_permutations(candidates_dir, gene_pdbs, output_dir)
