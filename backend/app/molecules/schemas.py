"""Pydantic schemas for molecule operations."""
from pydantic import BaseModel, Field, field_validator
from typing import Literal


class MoleculeConvertRequest(BaseModel):
    """Request to convert molecule between formats."""

    input_format: Literal["smiles", "pdb", "sdf", "mol2", "xyz"]
    output_format: Literal["smiles", "pdb", "sdf", "mol2", "xyz"]
    data: str = Field(..., description="Molecule data in input format")

    @field_validator("data")
    @classmethod
    def validate_data_not_empty(cls, v: str) -> str:
        """Validate that molecule data is not empty."""
        if not v.strip():
            raise ValueError("Molecule data cannot be empty")
        return v


class MoleculeInfo(BaseModel):
    """Molecular information."""

    formula: str
    molecular_weight: float
    num_atoms: int
    num_bonds: int
    num_heavy_atoms: int = 0
    num_rotatable_bonds: int = 0
    logp: float = 0.0


class MoleculeConvertResponse(BaseModel):
    """Response from molecule conversion."""

    output_format: str
    data: str
    molecule_info: MoleculeInfo


class MoleculeValidateRequest(BaseModel):
    """Request to validate molecule structure."""

    smiles: str = Field(..., description="SMILES string to validate")


class MoleculeValidateResponse(BaseModel):
    """Response from molecule validation."""

    valid: bool
    smiles: str
    canonical_smiles: str | None = None
    error_message: str | None = None


class MoleculePropertiesRequest(BaseModel):
    """Request to calculate molecular properties."""

    smiles: str = Field(..., description="SMILES string of molecule")


class MoleculePropertiesResponse(BaseModel):
    """Response with molecular properties."""

    smiles: str
    properties: MoleculeInfo
