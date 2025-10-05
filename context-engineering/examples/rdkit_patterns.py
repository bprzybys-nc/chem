"""
RDKit Integration Patterns

Real examples from chem platform showing RDKit integration patterns.
Based on backend/app/molecules/service.py (MoleculeService).
"""

from typing import Optional
from rdkit import Chem
from rdkit.Chem import Descriptors


def fast_failure_molecule_parsing(smiles: str) -> Chem.Mol:
    """
    Fast-failure pattern for molecule parsing.

    Pattern from MoleculeService.from_smiles():
    - Parse with RDKit
    - Check for None (invalid input)
    - Raise ValueError with actionable message
    - Include troubleshooting URL/guidance

    Args:
        smiles: SMILES string to parse

    Returns:
        RDKit molecule object

    Raises:
        ValueError: If SMILES is invalid with troubleshooting guidance
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            f"Invalid SMILES string: {smiles}. "
            f"Check syntax at: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"
        )
    return mol


def molecule_format_conversion_pattern(pdb_data: str) -> Chem.Mol:
    """
    Format conversion with troubleshooting guidance.

    Pattern from MoleculeService.from_pdb():
    - Parse with format-specific RDKit function
    - Check for None (invalid data)
    - Raise ValueError with multi-line troubleshooting
    - List specific checks for user to verify

    Args:
        pdb_data: PDB format string

    Returns:
        RDKit molecule object

    Raises:
        ValueError: With step-by-step troubleshooting
    """
    mol = Chem.MolFromPDBBlock(pdb_data)
    if mol is None:
        raise ValueError(
            "Invalid PDB data. "
            "Troubleshooting:\n"
            "  1. Verify PDB format is correct\n"
            "  2. Check for missing atom coordinates\n"
            "  3. Ensure file contains ATOM or HETATM records"
        )
    return mol


def property_calculation_pattern(mol: Chem.Mol) -> dict:
    """
    Property calculation with descriptive dict.

    Pattern from MoleculeService.calculate_properties():
    - Use RDKit Descriptors module
    - Return dict with descriptive keys
    - All values are simple JSON-serializable types
    - No abbreviations in keys (e.g., "molecular_weight" not "mw")

    Args:
        mol: RDKit molecule object

    Returns:
        Dict with molecular properties
    """
    return {
        "molecular_weight": Descriptors.MolWt(mol),
        "formula": Descriptors.rdMolDescriptors.CalcMolFormula(mol),
        "logp": Descriptors.MolLogP(mol),
        "num_atoms": mol.GetNumAtoms(),
        "num_bonds": mol.GetNumBonds(),
        "num_heavy_atoms": mol.GetNumHeavyAtoms(),
        "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
    }


def static_service_class_pattern():
    """
    Static service class pattern.

    Pattern from MoleculeService:
    - All methods are @staticmethod
    - No instance state (stateless operations)
    - Grouped by functionality (from_*, to_*, calculate_*)
    - Type hints on all parameters and returns
    - Comprehensive docstrings with Args/Returns/Raises

    Benefits:
    - No need to instantiate service
    - Clear separation of concerns
    - Easy to test (no mocking needed)
    - Follows functional programming principles
    """

    class MoleculeService:
        """Stateless service for molecule operations."""

        @staticmethod
        def from_smiles(smiles: str) -> Optional[Chem.Mol]:
            """Convert SMILES to molecule with fast failure."""
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {smiles}")
            return mol

        @staticmethod
        def to_pdb(mol: Chem.Mol, add_hydrogens: bool = True) -> str:
            """Convert molecule to PDB with optional hydrogens."""
            if add_hydrogens:
                mol = Chem.AddHs(mol)
            return Chem.MolToPDBBlock(mol)
