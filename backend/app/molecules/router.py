"""API router for molecule operations."""
from fastapi import APIRouter, HTTPException

from app.molecules import schemas, service

router = APIRouter()


@router.post("/convert", response_model=schemas.MoleculeConvertResponse)
async def convert_molecule(request: schemas.MoleculeConvertRequest):
    """
    Convert molecule between formats.

    Supports: SMILES, PDB, SDF, MOL2, XYZ

    Args:
        request: Conversion request with input/output formats and data

    Returns:
        Converted molecule data with properties

    Raises:
        HTTPException: If conversion fails
    """
    try:
        output_data, properties = service.MoleculeService.convert(
            request.input_format,
            request.output_format,
            request.data,
        )

        return schemas.MoleculeConvertResponse(
            output_format=request.output_format,
            data=output_data,
            molecule_info=schemas.MoleculeInfo(**properties),
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Conversion failed: {str(e)}. Please check input data and try again.",
        )


@router.post("/validate", response_model=schemas.MoleculeValidateResponse)
async def validate_molecule(request: schemas.MoleculeValidateRequest):
    """
    Validate SMILES string.

    Args:
        request: Validation request with SMILES string

    Returns:
        Validation result with canonical SMILES if valid

    Raises:
        HTTPException: If validation encounters unexpected error
    """
    try:
        mol = service.MoleculeService.from_smiles(request.smiles)
        canonical_smiles = service.MoleculeService.to_smiles(mol, canonical=True)

        return schemas.MoleculeValidateResponse(
            valid=True,
            smiles=request.smiles,
            canonical_smiles=canonical_smiles,
        )
    except ValueError as e:
        return schemas.MoleculeValidateResponse(
            valid=False,
            smiles=request.smiles,
            error_message=str(e),
        )
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Validation failed: {str(e)}",
        )


@router.post("/properties", response_model=schemas.MoleculePropertiesResponse)
async def calculate_properties(request: schemas.MoleculePropertiesRequest):
    """
    Calculate molecular properties from SMILES.

    Args:
        request: Request with SMILES string

    Returns:
        Molecular properties (formula, MW, atom counts, LogP, etc.)

    Raises:
        HTTPException: If calculation fails
    """
    try:
        mol = service.MoleculeService.from_smiles(request.smiles)
        properties = service.MoleculeService.calculate_properties(mol)

        return schemas.MoleculePropertiesResponse(
            smiles=request.smiles,
            properties=schemas.MoleculeInfo(**properties),
        )
    except ValueError as e:
        raise HTTPException(status_code=422, detail=str(e))
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Property calculation failed: {str(e)}",
        )
