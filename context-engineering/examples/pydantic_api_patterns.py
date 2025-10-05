"""
Pydantic API Schema Patterns

Real examples from chem platform showing Pydantic v2 schema design.
Based on backend/app/molecules/schemas.py.
"""

from pydantic import BaseModel, Field, field_validator


class MoleculeConvertRequest(BaseModel):
    """
    Request schema pattern with validation.

    Pattern from schemas.py:
    - Use Field() for descriptions and constraints
    - snake_case field names
    - Type hints for all fields
    - Optional validation with field_validator

    Example usage:
        request = MoleculeConvertRequest(
            input_format="smiles",
            output_format="pdb",
            data="c1ccccc1"
        )
    """

    input_format: str = Field(
        ...,
        description="Input molecule format (smiles, pdb, sdf, mol2, xyz)"
    )
    output_format: str = Field(
        ...,
        description="Output molecule format (smiles, pdb, sdf, mol2, xyz)"
    )
    data: str = Field(
        ...,
        description="Molecule data in input format"
    )

    @field_validator("data")
    @classmethod
    def validate_data(cls, v: str) -> str:
        """Fast-failure validation for required fields."""
        if not v or not v.strip():
            raise ValueError("Molecule data cannot be empty")
        return v


class MoleculeConvertResponse(BaseModel):
    """
    Response schema pattern with nested models.

    Pattern from schemas.py:
    - Descriptive field names
    - Optional fields for conditional data
    - Nested models for complex data (MoleculeInfo)
    - model_config for JSON examples
    """

    output_data: str = Field(
        ...,
        description="Converted molecule data"
    )
    molecule_info: "MoleculeInfo" = Field(
        ...,
        description="Molecular properties"
    )

    model_config = {
        "json_schema_extra": {
            "examples": [
                {
                    "output_data": "ATOM    1  C   ...",
                    "molecule_info": {
                        "molecular_weight": 78.11,
                        "formula": "C6H6"
                    }
                }
            ]
        }
    }


class MoleculeInfo(BaseModel):
    """
    Nested model pattern for complex data.

    Pattern from schemas.py:
    - Simple, focused responsibility
    - All fields have descriptions
    - Optional fields use Optional[type] or None default
    - Used within other response schemas
    """

    molecular_weight: float = Field(..., description="Molecular weight in g/mol")
    formula: str = Field(..., description="Molecular formula")
    logp: float = Field(..., description="Partition coefficient (LogP)")
    num_atoms: int = Field(..., description="Total number of atoms")
    num_bonds: int = Field(..., description="Total number of bonds")
    num_heavy_atoms: int = Field(..., description="Number of non-hydrogen atoms")
    num_rotatable_bonds: int = Field(..., description="Number of rotatable bonds")


class MoleculeValidateResponse(BaseModel):
    """
    Validation response pattern with conditional fields.

    Pattern from schemas.py:
    - Boolean success flag
    - Optional success data (canonical_smiles)
    - Optional error data (error_message)
    - Client can check 'valid' field to determine which fields are present
    """

    valid: bool = Field(..., description="Whether SMILES is valid")
    smiles: str = Field(..., description="Original SMILES string")
    canonical_smiles: str | None = Field(
        None,
        description="Canonical SMILES form (if valid)"
    )
    error_message: str | None = Field(
        None,
        description="Error message (if invalid)"
    )


# Pattern: Separation of Request/Response schemas
# Each API endpoint has dedicated request and response schemas
# Benefits:
# - Clear API contract
# - Type-safe validation
# - Auto-generated OpenAPI docs
# - Easy client generation
