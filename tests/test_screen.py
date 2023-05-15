#!/usr/bin/env python3

import pytest
import tempfile
from pathlib import Path

import screen
import util
import rcsb_pdb
import chem


@pytest.fixture
def data_dir(): return Path(__file__).parent / "data"

def test_dock_fabp4_with_bms309403_from_smiles(data_dir):
    chem.compounds()
    receptor_dir = util.mkdir(screen.data_dir() / 'receptors/3RZY.pdb')
    output_dir = util.mkdir(data_dir / 'outputs')
    def dock_to_fabp4(compound, output_stem): return screen.dock(compound, receptor_dir, output_dir / output_stem)
    breakpoint()
    smiles_docking = dock_to_fabp4(util.mkdir(screen.data_dir() / 'compounds') / 'BMS.sdf',
                                     'bms_from_smiles_with_3rzy')
    rcsb_docking = dock_to_fabp4(screen.pdb_partition(rcsb_pdb.write_pdb("2NNQ"), "T4B"),
                                 'bms_from_rcsb_with_3rzy')
    assert util.util.files_are_equivalent(smiles_docking["sdf"], rcsb_docking["sdf"])


if __name__=="__main__":
    test_dock_fabp4_with_bms309403_from_smiles()
