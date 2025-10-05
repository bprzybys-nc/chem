"""Unit tests for molecule operations."""
import pytest


class TestHealthEndpoint:
    """Test health check endpoint."""

    def test_health_check(self, client):
        """Test health endpoint returns healthy status."""
        response = client.get("/health")
        assert response.status_code == 200
        data = response.json()
        assert data["status"] == "healthy"
        assert data["service"] == "chem-backend"


class TestMoleculeValidation:
    """Test molecule validation endpoint."""

    def test_validate_valid_smiles(self, client, benzene_smiles):
        """Test validation with valid SMILES."""
        response = client.post(
            "/api/molecules/validate",
            json={"smiles": benzene_smiles},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
        assert data["smiles"] == benzene_smiles
        assert data["canonical_smiles"] is not None

    def test_validate_invalid_smiles(self, client, invalid_smiles):
        """Test validation with invalid SMILES."""
        response = client.post(
            "/api/molecules/validate",
            json={"smiles": invalid_smiles},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is False
        assert data["error_message"] is not None
        assert "Invalid SMILES" in data["error_message"]

    def test_validate_empty_smiles(self, client):
        """Test validation with empty SMILES."""
        response = client.post(
            "/api/molecules/validate",
            json={"smiles": ""},
        )
        # RDKit accepts empty string but returns empty molecule
        assert response.status_code == 200
        data = response.json()
        # Empty SMILES is technically valid for RDKit (creates empty molecule)
        assert data["valid"] is True


class TestMoleculeConversion:
    """Test molecule format conversion endpoint."""

    def test_convert_smiles_to_pdb(self, client, benzene_smiles):
        """Test SMILES to PDB conversion."""
        response = client.post(
            "/api/molecules/convert",
            json={
                "input_format": "smiles",
                "output_format": "pdb",
                "data": benzene_smiles,
            },
        )
        assert response.status_code == 200
        data = response.json()
        assert data["output_format"] == "pdb"
        assert "HETATM" in data["data"] or "ATOM" in data["data"]
        assert data["molecule_info"]["formula"] == "C6H6"
        # PDB conversion adds hydrogens (to_pdb has add_hydrogens=True)
        assert data["molecule_info"]["num_atoms"] == 12

    def test_convert_pdb_to_smiles(self, client, benzene_pdb):
        """Test PDB to SMILES conversion."""
        response = client.post(
            "/api/molecules/convert",
            json={
                "input_format": "pdb",
                "output_format": "smiles",
                "data": benzene_pdb,
            },
        )
        assert response.status_code == 200
        data = response.json()
        assert data["output_format"] == "smiles"
        assert len(data["data"]) > 0
        # Should contain benzene ring
        assert "c" in data["data"].lower()

    def test_convert_invalid_smiles(self, client, invalid_smiles):
        """Test conversion with invalid SMILES returns error."""
        response = client.post(
            "/api/molecules/convert",
            json={
                "input_format": "smiles",
                "output_format": "pdb",
                "data": invalid_smiles,
            },
        )
        assert response.status_code == 422
        assert "Invalid SMILES" in response.json()["detail"]

    def test_convert_empty_data(self, client):
        """Test conversion with empty data returns error."""
        response = client.post(
            "/api/molecules/convert",
            json={
                "input_format": "smiles",
                "output_format": "pdb",
                "data": "   ",  # Only whitespace
            },
        )
        assert response.status_code == 422
        # Pydantic v2 returns error list
        detail = response.json()["detail"]
        if isinstance(detail, list):
            assert "cannot be empty" in str(detail)
        else:
            assert "cannot be empty" in detail


class TestMoleculeProperties:
    """Test molecular property calculation endpoint."""

    def test_calculate_properties(self, client, benzene_smiles):
        """Test property calculation for benzene."""
        response = client.post(
            "/api/molecules/properties",
            json={"smiles": benzene_smiles},
        )
        assert response.status_code == 200
        data = response.json()
        assert data["smiles"] == benzene_smiles

        props = data["properties"]
        assert props["formula"] == "C6H6"
        assert 78.0 < props["molecular_weight"] < 79.0  # ~78.11
        # RDKit SMILES parsing doesn't add hydrogens by default
        # Aromatic benzene has 6 carbons only
        assert props["num_atoms"] == 6  # 6 carbons (implicit hydrogens)
        assert props["num_bonds"] == 6  # Aromatic bonds
        assert props["num_heavy_atoms"] == 6  # Only carbons

    def test_calculate_properties_invalid_smiles(self, client, invalid_smiles):
        """Test property calculation with invalid SMILES."""
        response = client.post(
            "/api/molecules/properties",
            json={"smiles": invalid_smiles},
        )
        assert response.status_code == 422
        assert "Invalid SMILES" in response.json()["detail"]


class TestAPIDocumentation:
    """Test OpenAPI documentation."""

    def test_openapi_docs_accessible(self, client):
        """Test that OpenAPI documentation is accessible."""
        response = client.get("/docs")
        assert response.status_code == 200

    def test_openapi_json_accessible(self, client):
        """Test that OpenAPI JSON schema is accessible."""
        response = client.get("/openapi.json")
        assert response.status_code == 200
        data = response.json()
        assert "openapi" in data
        assert "paths" in data
        assert "/api/molecules/convert" in data["paths"]
        assert "/api/molecules/validate" in data["paths"]
        assert "/api/molecules/properties" in data["paths"]
