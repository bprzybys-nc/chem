"""Service for molecule operations using RDKit."""
from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors


class MoleculeService:
    """Service for molecule operations using RDKit."""

    @staticmethod
    def from_smiles(smiles: str) -> Optional[Chem.Mol]:
        """
        Convert SMILES to RDKit molecule.

        Fast-failure validation - raises ValueError with actionable message.

        Args:
            smiles: SMILES string

        Returns:
            RDKit molecule object

        Raises:
            ValueError: If SMILES string is invalid
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(
                f"Invalid SMILES string: {smiles}. "
                f"Check syntax at: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"
            )
        return mol

    @staticmethod
    def from_pdb(pdb_data: str) -> Optional[Chem.Mol]:
        """
        Convert PDB data to RDKit molecule.

        Args:
            pdb_data: PDB format data

        Returns:
            RDKit molecule object

        Raises:
            ValueError: If PDB data is invalid
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

    @staticmethod
    def from_mol_block(mol_data: str, format_type: str) -> Optional[Chem.Mol]:
        """
        Convert MOL/SDF data to RDKit molecule.

        Args:
            mol_data: MOL/SDF format data
            format_type: Either 'sdf' or 'mol2'

        Returns:
            RDKit molecule object

        Raises:
            ValueError: If molecule data is invalid
        """
        mol = Chem.MolFromMolBlock(mol_data)
        if mol is None:
            raise ValueError(
                f"Invalid {format_type.upper()} data. "
                "Troubleshooting:\n"
                "  1. Verify file format is correct\n"
                "  2. Check for valid atom/bond records\n"
                "  3. Ensure coordinate data is present"
            )
        return mol

    @staticmethod
    def from_xyz(xyz_data: str) -> Optional[Chem.Mol]:
        """
        Convert XYZ data to RDKit molecule.

        Args:
            xyz_data: XYZ format data

        Returns:
            RDKit molecule object

        Raises:
            ValueError: If XYZ data is invalid
        """
        mol = Chem.MolFromXYZBlock(xyz_data)
        if mol is None:
            raise ValueError(
                "Invalid XYZ data. "
                "Troubleshooting:\n"
                "  1. Verify XYZ format (atom count, comment, coordinates)\n"
                "  2. Check coordinate values are valid numbers\n"
                "  3. Ensure element symbols are correct"
            )
        return mol

    @staticmethod
    def to_smiles(mol: Chem.Mol, canonical: bool = True) -> str:
        """Convert RDKit molecule to SMILES."""
        if canonical:
            return Chem.MolToSmiles(mol)
        return Chem.MolToSmiles(mol, canonical=False)

    @staticmethod
    def to_pdb(mol: Chem.Mol, add_hydrogens: bool = True) -> str:
        """
        Convert RDKit molecule to PDB format.

        Args:
            mol: RDKit molecule object
            add_hydrogens: Whether to add hydrogen atoms

        Returns:
            PDB format string
        """
        if add_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

        return Chem.MolToPDBBlock(mol)

    @staticmethod
    def to_mol_block(mol: Chem.Mol) -> str:
        """Convert RDKit molecule to MOL block (SDF format)."""
        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

        return Chem.MolToMolBlock(mol)

    @staticmethod
    def to_xyz(mol: Chem.Mol) -> str:
        """Convert RDKit molecule to XYZ format."""
        # Generate 3D coordinates if not present
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

        return Chem.MolToXYZBlock(mol)

    @staticmethod
    def calculate_properties(mol: Chem.Mol) -> dict:
        """
        Calculate molecular properties.

        Args:
            mol: RDKit molecule object

        Returns:
            Dictionary of molecular properties
        """
        return {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": Descriptors.MolWt(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "logp": Descriptors.MolLogP(mol),
        }

    @staticmethod
    def convert(input_format: str, output_format: str, data: str) -> tuple[str, dict]:
        """
        Convert molecule between formats.

        Args:
            input_format: Input format (smiles, pdb, sdf, mol2, xyz)
            output_format: Output format (smiles, pdb, sdf, mol2, xyz)
            data: Molecule data in input format

        Returns:
            Tuple of (converted_data, properties_dict)

        Raises:
            ValueError: If conversion fails
        """
        # Parse input format
        if input_format == "smiles":
            mol = MoleculeService.from_smiles(data)
        elif input_format == "pdb":
            mol = MoleculeService.from_pdb(data)
        elif input_format in ("sdf", "mol2"):
            mol = MoleculeService.from_mol_block(data, input_format)
        elif input_format == "xyz":
            mol = MoleculeService.from_xyz(data)
        else:
            raise ValueError(f"Unsupported input format: {input_format}")

        # Convert to output format
        if output_format == "smiles":
            output_data = MoleculeService.to_smiles(mol)
        elif output_format == "pdb":
            output_data = MoleculeService.to_pdb(mol)
            # PDB conversion adds hydrogens, update mol for property calculation
            mol = Chem.AddHs(mol)
        elif output_format in ("sdf", "mol2"):
            output_data = MoleculeService.to_mol_block(mol)
        elif output_format == "xyz":
            output_data = MoleculeService.to_xyz(mol)
        else:
            raise ValueError(f"Unsupported output format: {output_format}")

        # Calculate properties after conversion (reflects any modifications like H addition)
        properties = MoleculeService.calculate_properties(mol)

        return output_data, properties
