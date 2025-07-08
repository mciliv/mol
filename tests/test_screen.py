#!/usr/bin/env python3

import tempfile
from pathlib import Path
import unittest
import filecmp

import pytest

import screen
import rcsb_pdb
import chem


def local_data_dir(file_path):
    """Create a local data directory based on the test file location."""
    test_dir = Path(file_path).parent
    data_dir = test_dir / "test_data"
    data_dir.mkdir(parents=True, exist_ok=True)
    return data_dir


def mkdirs(path):
    """Create directories ensuring the path exists."""
    path = Path(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def files_are_equivalent(file1, file2):
    """Compare two files to see if they are equivalent."""
    try:
        return filecmp.cmp(file1, file2, shallow=False)
    except (OSError, FileNotFoundError):
        return False


# @pytest.fixture
def local_data_dir_fixture():
    return local_data_dir(__file__)


def test_dock_fabp4_with_bms309403_from_smiles(local_data_dir=None):
    if local_data_dir is None:
        local_data_dir = local_data_dir_fixture()
    
    receptor = mkdirs(screen.local_data_dir / "receptors") / "3RZY.pdb"
    output_dir = mkdirs(local_data_dir / "outputs")

    def dock_to_fabp4(compound, output_stem):
        return screen.dock(receptor, compound, output_dir / output_stem)

    smiles_docking = dock_to_fabp4(
        chem.sdf(
            "BMS",
            "[O-][S](=O)(=O)c1cccc2cccc(Nc3ccccc3)c12",
            directory_path=local_data_dir / "compounds",
        ),
        "bms_from_smiles_with_3rzy",
    )
    rcsb_docking = dock_to_fabp4(
        screen.pdb_partition(rcsb_pdb.write_pdb("2NNQ"), "T4B"),
        "bms_from_rcsb_with_3rzy",
    )
    assert files_are_equivalent(smiles_docking["sdf"], rcsb_docking["sdf"])


if __name__ == "__main__":
    test_dock_fabp4_with_bms309403_from_smiles(local_data_dir_fixture())
