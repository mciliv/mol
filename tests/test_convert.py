import pytest
from pathlib import Path
from convert.sdf import sdf

@pytest.fixture
def cleanup_sdf():
    """Fixture to clean up any generated SDF files after the test."""
    yield
    for sdf_file in Path("test_sdf_output").glob("*.sdf"):
        sdf_file.unlink()
    Path("test_sdf_output").rmdir()

@pytest.mark.parametrize("smiles", ["CCO", "O", "C1=CC=CC=C1"])
def test_sdf_file_creation(smiles, cleanup_sdf):
    directory = "test_sdf_output"
    Path(directory).mkdir(exist_ok=True)

    filepath = sdf(smiles, directory, overwrite=True)
    
    assert Path(filepath).exists(), f"SDF file {filepath} was not created."
