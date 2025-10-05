# Project Overview

## Purpose
Chemistry 3D visualization platform - web application for visualizing and analyzing molecular structures in 3D.

## Tech Stack

### Backend
- **Framework**: FastAPI (async web framework)
- **Chemistry**: RDKit (molecule processing, format conversion, property calculation)
- **Validation**: Pydantic v2 (data validation and type safety)
- **Server**: Uvicorn
- **Package Manager**: UV (uv package management)
- **Python**: >=3.11

### Frontend  
- **Framework**: Next.js 15.5.4 (React 19.1.0)
- **Visualization**: 3Dmol.js (3D molecule rendering)
- **Styling**: Tailwind CSS v4
- **API Client**: @hey-api/client-fetch (generated from OpenAPI)
- **TypeScript**: Full type safety

### Infrastructure
- **Containerization**: Docker + Docker Compose
- **Backend Port**: 8000
- **Frontend Port**: 3000

## Project Structure
```
chem/
├── backend/               # FastAPI + RDKit backend (477 LOC)
│   ├── app/
│   │   ├── main.py       # FastAPI entry point (48 lines)
│   │   ├── config.py     # Configuration (26 lines)
│   │   └── molecules/    # Molecule operations module
│   │       ├── router.py     # API endpoints (110 lines)
│   │       ├── schemas.py    # Pydantic models (67 lines, 7 classes)
│   │       └── service.py    # RDKit operations (226 lines, 10 methods)
│   ├── tests/
│   │   ├── conftest.py
│   │   └── unit/
│   │       └── test_molecules.py  # 12 tests passing
│   └── pyproject.toml
│
├── frontend/              # Next.js + 3Dmol.js frontend (221 LOC)
│   ├── app/
│   │   ├── page.tsx              # Landing page (35 lines)
│   │   └── visualize/
│   │       └── page.tsx          # Visualization page (114 lines)
│   ├── components/
│   │   └── MoleculeViewer.tsx    # 3D viewer (72 lines)
│   ├── src/api/                  # Generated API client
│   ├── types/
│   └── package.json
│
├── context-engineering/   # PRP methodology framework
│   ├── PRPs/              # 1 validated PRP (PRP-0)
│   ├── examples/          # Code pattern examples
│   └── templates/         # INITIAL templates
│
└── docker-compose.yml

## Current Implementation Status

### Phase 1: Backend Foundation ✅ COMPLETE
**Implemented Components:**
- `MoleculeService` class (service.py) - 10 static methods
  - Format parsers: `from_smiles()`, `from_pdb()`, `from_mol_block()`, `from_xyz()`
  - Format converters: `to_smiles()`, `to_pdb()`, `to_mol_block()`, `to_xyz()`
  - Analysis: `calculate_properties()`
  - High-level: `convert()` wrapper
  
- Pydantic Schemas (schemas.py) - 7 classes
  - `MoleculeConvertRequest`, `MoleculeConvertResponse`
  - `MoleculeValidateRequest`, `MoleculeValidateResponse`
  - `MoleculePropertiesRequest`, `MoleculePropertiesResponse`
  - `MoleculeInfo` (nested model)

- FastAPI Router (router.py) - 3 endpoints
  - `POST /api/molecules/validate` - Validate SMILES
  - `POST /api/molecules/convert` - Convert formats
  - `POST /api/molecules/properties` - Calculate properties
  - Plus: `GET /health` in main.py

- React Components (frontend)
  - `MoleculeViewer` - 3D visualization with 3Dmol.js
  - Dynamic import pattern (SSR disabled for WebGL)

**Testing:**
- 12/12 unit tests passing
- Real functionality testing (no mocks)
- FastAPI TestClient integration

**Quality Metrics:**
- All files under 500-line guideline ✅
- All functions under 50-line guideline ✅
- 100% CLAUDE.md compliance ✅
- Fast-failure error handling ✅

## Features Implemented

✅ Molecule format conversion (SMILES ↔ PDB, SDF, MOL2, XYZ)
✅ Molecule validation (SMILES with canonical form)
✅ Property calculation (MW, formula, LogP, atoms, bonds, etc.)
✅ RESTful API with auto-generated OpenAPI docs
✅ Health check endpoint
✅ 3D molecule visualization (basic)
✅ Type-safe API client generation
✅ Docker containerization
✅ Fast-failure error handling with troubleshooting

## Next Steps (Future Phases)

**Phase 2 Potential Features:**
- Enhanced 3D visualization controls
- Multiple molecule display styles
- File upload for molecule formats
- Molecule comparison view
- Property visualization

## Key Design Patterns

- **Service Layer Pattern**: Stateless service with static methods
- **Fast-Failure Pattern**: Immediate error raising with guidance
- **Schema/DTO Pattern**: Separate request/response models
- **Dynamic Import Pattern**: Client-side only library loading
- **Generated Client Pattern**: OpenAPI → TypeScript client