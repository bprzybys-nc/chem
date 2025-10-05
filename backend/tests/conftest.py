"""Pytest fixtures and configuration."""
import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture
def client():
    """FastAPI test client fixture."""
    return TestClient(app)


@pytest.fixture
def benzene_smiles():
    """Benzene SMILES string fixture."""
    return "c1ccccc1"


@pytest.fixture
def benzene_pdb():
    """Benzene PDB data fixture."""
    return """COMPND    benzene
HETATM    1  C   UNL     1       1.200   0.000   0.000  1.00  0.00           C
HETATM    2  C   UNL     1       0.600   1.039   0.000  1.00  0.00           C
HETATM    3  C   UNL     1      -0.600   1.039   0.000  1.00  0.00           C
HETATM    4  C   UNL     1      -1.200   0.000   0.000  1.00  0.00           C
HETATM    5  C   UNL     1      -0.600  -1.039   0.000  1.00  0.00           C
HETATM    6  C   UNL     1       0.600  -1.039   0.000  1.00  0.00           C
HETATM    7  H   UNL     1       2.130   0.000   0.000  1.00  0.00           H
HETATM    8  H   UNL     1       1.065   1.845   0.000  1.00  0.00           H
HETATM    9  H   UNL     1      -1.065   1.845   0.000  1.00  0.00           H
HETATM   10  H   UNL     1      -2.130   0.000   0.000  1.00  0.00           H
HETATM   11  H   UNL     1      -1.065  -1.845   0.000  1.00  0.00           H
HETATM   12  H   UNL     1       1.065  -1.845   0.000  1.00  0.00           H
CONECT    1    2    6    7
CONECT    2    1    3    8
CONECT    3    2    4    9
CONECT    4    3    5   10
CONECT    5    4    6   11
CONECT    6    1    5   12
END
"""


@pytest.fixture
def invalid_smiles():
    """Invalid SMILES string fixture."""
    return "invalid_smiles_xyz"
