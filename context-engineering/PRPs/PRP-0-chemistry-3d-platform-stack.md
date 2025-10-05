---
name: "Chemistry 3D Visualization Platform - Tech Stack"
description: "Complete tech stack blueprint for FastAPI backend + Next.js frontend with 3D molecule visualization using RDKit and 3Dmol.js"
prp_id: "PRP-0"
status: "new"
created_date: "2025-10-05T22:45:00Z"
last_updated: "2025-10-05T23:00:00Z"
updated_by: "peer-review-command"
context_sync:
  ce_updated: false
  serena_updated: false
version: 1.1
---

# Chemistry 3D Visualization Platform - Tech Stack

## Context

Building a modern web platform for interactive 3D visualization of chemical molecules. The system will:
- Accept multiple chemistry data formats (SMILES, PDB, SDF, MOL2, XYZ)
- Process and convert molecules using industry-standard chemistry libraries
- Render interactive 3D visualizations in the browser
- Provide RESTful API for molecule operations
- Support both desktop and mobile devices

**Key Challenge:** Integrate specialized chemistry processing (RDKit) with modern web technologies while maintaining performance and type safety.

**Target Users:** Chemists, researchers, educators, students working with molecular structures.

## Requirements

### Functional Requirements
1. **Molecule Input**
   - Accept SMILES strings
   - Upload PDB, SDF, MOL2, XYZ files
   - Validate molecule structure

2. **Molecule Processing**
   - Convert between formats
   - Calculate molecular properties (MW, formula, etc.)
   - Perform substructure searches
   - Add/remove hydrogens

3. **3D Visualization**
   - Interactive 3D rendering (rotate, zoom, pan)
   - Multiple representation styles (ball-stick, spacefill, cartoon)
   - Color by atom type, element, property
   - Support molecules with 500+ atoms

4. **API Layer**
   - RESTful endpoints with OpenAPI documentation
   - Type-safe client generation
   - JSON request/response (Chemical JSON format)
   - Error handling with actionable messages

### Non-Functional Requirements
1. **Performance**
   - API response < 200ms for simple molecules
   - 60 FPS rendering for structures with 500+ atoms
   - Progressive loading for large molecules

2. **Developer Experience**
   - Full TypeScript type safety
   - Hot reload in development
   - Automated API client generation
   - Comprehensive testing

3. **Deployment**
   - Containerized with Docker
   - Environment-based configuration
   - Health checks for all services
   - Production-ready logging

## Tech Stack

### Frontend Stack

**Framework & Language:**
- **Next.js 14+** - React framework with App Router
- **TypeScript** - Type-safe JavaScript
- **React 18+** - UI library

**3D Visualization:**
- **3Dmol.js** - WebGL-based molecular visualization
  - Most popular option (847 GitHub stars)
  - Supports PDB, SDF, MOL2, XYZ, CIF, MMTF formats
  - BSD open-source license
  - Easy embedding with minimal code
- **molecule-3d-for-react** - React wrapper for 3Dmol.js
  - Data-bound React component
  - Maintained by Autodesk

**API Integration:**
- **@hey-api/openapi-ts** - Type-safe API client generation from OpenAPI spec
- **Axios** or **Fetch API** - HTTP client

**Styling:**
- **Tailwind CSS** - Utility-first CSS framework
- **shadcn/ui** - Component library (optional)

### Backend Stack

**Framework & Language:**
- **FastAPI** - Modern async Python web framework
  - Auto-generated OpenAPI documentation
  - Pydantic v2 integration
  - High performance (comparable to Node.js)
  - Excellent for data APIs

**Chemistry Processing:**
- **RDKit** - Open-source cheminformatics toolkit
  - SMILES parsing and generation
  - Format conversion (PDB, SDF, MOL2, XYZ)
  - Molecular property calculations
  - Substructure search
  - Fingerprints and similarity
  - Pandas integration

**Data Validation:**
- **Pydantic v2** - Data validation and serialization
  - Type hints for validation
  - JSON schema generation
  - Automatic OpenAPI integration

**Project Structure:**
- Module-based organization (not file-type)
- Separation of concerns (router, schema, service, utils)

### Data Formats

**Input Formats:**
- **SMILES** - Simplified Molecular-Input Line-entry System (string representation)
- **PDB** - Protein Data Bank (3D structures, most common for proteins)
- **SDF** - Structure-Data File (supports metadata and properties)
- **MOL2** - Tripos format (includes partial charges)
- **XYZ** - Simple coordinates format (basic, no connectivity)

**API Response Format:**
- **Chemical JSON** - Compact JSON format for molecular data
  - Arrays for atoms, bonds, coordinates
  - 3D positions offset by 3N (N = atom index)
  - Efficient for web transfer

**Format Conversion Flow:**
```
User Input (any format) → RDKit → Chemical JSON → Frontend → 3Dmol.js
```

### Deployment Stack

**Containerization:**
- **Docker** - Container platform
- **Docker Compose** - Multi-service orchestration

**Deployment Platforms:**
- **Railway.app** - Recommended for managed deployment
- **Render** - Alternative managed platform
- Both support Docker containers natively

**Architecture:**
- Monorepo with backend/ and frontend/ directories
- Separate Dockerfiles for each service
- docker-compose.yml for local development

## Project Structure

```
chem/
├── backend/                          # FastAPI application
│   ├── app/
│   │   ├── __init__.py
│   │   ├── main.py                  # FastAPI app entry point
│   │   ├── config.py                # Settings with Pydantic BaseSettings
│   │   │
│   │   ├── molecules/               # Module: molecule operations
│   │   │   ├── __init__.py
│   │   │   ├── router.py            # API endpoints (/api/molecules/*)
│   │   │   ├── schemas.py           # Pydantic models (request/response)
│   │   │   ├── service.py           # Business logic with RDKit
│   │   │   └── utils.py             # Helper functions
│   │   │
│   │   └── visualization/           # Module: visualization data prep
│   │       ├── __init__.py
│   │       ├── router.py            # API endpoints (/api/visualization/*)
│   │       ├── schemas.py           # Pydantic models
│   │       └── service.py           # Chemical JSON generation
│   │
│   ├── tests/
│   │   ├── unit/                    # Fast, isolated tests
│   │   └── integration/             # API integration tests
│   │
│   ├── pyproject.toml               # UV package management
│   ├── Dockerfile                   # Backend container
│   └── .env.example                 # Environment template
│
├── frontend/                         # Next.js application
│   ├── app/                         # App Router
│   │   ├── layout.tsx               # Root layout
│   │   ├── page.tsx                 # Home page
│   │   └── visualize/
│   │       └── page.tsx             # Visualization page
│   │
│   ├── components/
│   │   ├── MoleculeViewer.tsx       # 3Dmol.js wrapper (client-only)
│   │   ├── MoleculeUpload.tsx       # File upload component
│   │   └── MoleculeProperties.tsx   # Property display
│   │
│   ├── lib/
│   │   └── api/                     # Generated API client
│   │       └── client.ts            # Type-safe API client
│   │
│   ├── types/
│   │   └── molecule.ts              # TypeScript types
│   │
│   ├── package.json
│   ├── tsconfig.json
│   ├── next.config.js
│   ├── Dockerfile                   # Frontend container
│   └── .env.example
│
├── docker-compose.yml               # Multi-service orchestration
├── .gitignore
├── README.md
└── context-engineering/             # Context Engineering framework
    ├── PRPs/
    │   └── chemistry-3d-platform-stack.md  # This file
    ├── examples/
    └── templates/
```

## API Design Patterns

### Endpoint Structure

Following REST best practices with plural nouns and HTTP verbs:

**Molecule Operations:**
```
POST   /api/molecules/convert        # Convert between formats
POST   /api/molecules/validate       # Validate molecule structure
POST   /api/molecules/properties     # Calculate molecular properties
POST   /api/molecules/search         # Substructure search
```

**Visualization Data:**
```
POST   /api/visualization/prepare    # Prepare molecule for 3D rendering
GET    /api/visualization/formats    # List supported formats
```

### Example API Requests/Responses

**Convert SMILES to PDB:**
```bash
curl -X POST http://localhost:8000/api/molecules/convert \
  -H "Content-Type: application/json" \
  -d '{
    "input_format": "smiles",
    "output_format": "pdb",
    "data": "c1ccccc1"
  }'
```

Response:
```json
{
  "output_format": "pdb",
  "data": "COMPND    benzene\nHETATM...",
  "molecule_info": {
    "formula": "C6H6",
    "molecular_weight": 78.11,
    "num_atoms": 12,
    "num_bonds": 12
  }
}
```

**Prepare for Visualization:**
```bash
curl -X POST http://localhost:8000/api/visualization/prepare \
  -H "Content-Type: application/json" \
  -d '{
    "smiles": "c1ccccc1",
    "style": "stick"
  }'
```

Response (Chemical JSON):
```json
{
  "chemical_json": {
    "atoms": {
      "elements": [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1],
      "coords_3d": [1.2, 0.0, 0.0, 0.6, 1.04, 0.0, ...],
      "formal_charges": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    },
    "bonds": {
      "connections": [[0, 1], [1, 2], [2, 3], [3, 4], [4, 5], [5, 0], ...],
      "order": [2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 1, 1]
    }
  },
  "pdb_data": "COMPND    benzene\nHETATM...",
  "properties": {
    "formula": "C6H6",
    "molecular_weight": 78.11
  }
}
```

### Pydantic Schema Examples

**molecules/schemas.py:**
```python
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
        if not v.strip():
            raise ValueError("Molecule data cannot be empty")
        return v

class MoleculeInfo(BaseModel):
    """Molecular information."""
    formula: str
    molecular_weight: float
    num_atoms: int
    num_bonds: int

class MoleculeConvertResponse(BaseModel):
    """Response from molecule conversion."""
    output_format: str
    data: str
    molecule_info: MoleculeInfo
```

**visualization/schemas.py:**
```python
from pydantic import BaseModel, model_validator
from typing import List, Optional, Union, Literal

class ChemicalJSON(BaseModel):
    """Chemical JSON format for molecular data."""
    atoms: dict[str, List[Union[float, int]]]
    bonds: dict[str, Union[List[List[int]], List[int]]]

class VisualizationPrepareRequest(BaseModel):
    """Request to prepare molecule for visualization."""
    smiles: Optional[str] = None
    pdb_data: Optional[str] = None
    style: Literal["stick", "sphere", "cartoon"] = "stick"

    @model_validator(mode='after')
    def check_at_least_one_input(self) -> 'VisualizationPrepareRequest':
        if not self.smiles and not self.pdb_data:
            raise ValueError("Either 'smiles' or 'pdb_data' must be provided")
        return self

class VisualizationPrepareResponse(BaseModel):
    """Response with visualization-ready data."""
    chemical_json: ChemicalJSON
    pdb_data: str
    properties: MoleculeInfo
```

### RDKit Service Examples

**molecules/service.py:**
```python
from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from typing import Optional

class MoleculeService:
    """Service for molecule operations using RDKit."""

    @staticmethod
    def from_smiles(smiles: str) -> Optional[Chem.Mol]:
        """
        Convert SMILES to RDKit molecule.

        Fast-failure validation - raises ValueError with actionable message.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(
                f"Invalid SMILES string: {smiles}. "
                f"Check syntax at: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html"
            )
        return mol

    @staticmethod
    def to_pdb(mol: Chem.Mol, add_hydrogens: bool = True) -> str:
        """Convert RDKit molecule to PDB format."""
        if add_hydrogens:
            mol = Chem.AddHs(mol)

        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        return Chem.MolToPDBBlock(mol)

    @staticmethod
    def calculate_properties(mol: Chem.Mol) -> dict:
        """Calculate molecular properties."""
        return {
            "formula": Chem.rdMolDescriptors.CalcMolFormula(mol),
            "molecular_weight": Descriptors.MolWt(mol),
            "num_atoms": mol.GetNumAtoms(),
            "num_bonds": mol.GetNumBonds(),
            "num_heavy_atoms": mol.GetNumHeavyAtoms(),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "logp": Descriptors.MolLogP(mol),
        }
```

## Implementation Approach

### Phase 1: Backend Foundation (Week 1)

**Setup:**
1. Initialize FastAPI project with UV
2. Install RDKit (use conda or Docker for consistent environment)
3. Configure module-based structure

**Implementation:**
1. Create `main.py` with FastAPI app and CORS middleware
2. Add health check endpoint at `/health`
3. Implement `molecules` module:
   - Router with `/convert`, `/validate`, `/properties` endpoints
   - Pydantic schemas for request/response
   - Service layer with RDKit operations
   - Error handling with fast-failure pattern
4. Add comprehensive unit tests with fixtures
5. Generate OpenAPI documentation

**Key Files:**
- `backend/app/main.py` - FastAPI entry point with CORS and health check
- `backend/app/molecules/router.py` - API endpoints
- `backend/app/molecules/schemas.py` - Pydantic models
- `backend/app/molecules/service.py` - RDKit operations
- `backend/tests/conftest.py` - Pytest fixtures and configuration
- `backend/tests/unit/test_molecules.py` - Unit tests

### Phase 2: Frontend Foundation (Week 2)

**Setup:**
1. Initialize Next.js project with TypeScript
2. Install 3Dmol.js and molecule-3d-for-react
3. Configure for client-side WebGL rendering

**Implementation:**
1. Create MoleculeViewer component:
   - Use dynamic import for client-side only
   - Wrap 3Dmol.js with React component
   - Implement interactive controls
2. Create basic UI for molecule input
3. Test 3D rendering with hardcoded molecules

**Key Files:**
- `frontend/components/MoleculeViewer.tsx` - 3Dmol.js wrapper
- `frontend/app/visualize/page.tsx` - Visualization page
- `frontend/lib/utils.ts` - Helper functions

**WebGL Optimization:**
```typescript
// MoleculeViewer.tsx - Client-side only component
'use client'

import dynamic from 'next/dynamic'
import { useEffect, useRef } from 'react'

interface MoleculeViewerProps {
  pdbData: string
}

// Dynamic import to avoid SSR
const MoleculeViewer: React.FC<MoleculeViewerProps> = ({ pdbData }) => {
  const viewerRef = useRef<HTMLDivElement>(null)
  const $3Dmol = useRef<any>(null)

  useEffect(() => {
    // Load 3Dmol.js only on client
    if (typeof window !== 'undefined' && !$3Dmol.current) {
      import('3dmol/build/3Dmol.js').then((module) => {
        $3Dmol.current = module.default
        initializeViewer()
      })
    }
  }, []) // Initial load only

  useEffect(() => {
    // Update viewer when pdbData changes
    if ($3Dmol.current && pdbData) {
      initializeViewer()
    }
  }, [pdbData])

  const initializeViewer = () => {
    if (!viewerRef.current || !$3Dmol.current) return

    const viewer = $3Dmol.current.createViewer(viewerRef.current, {
      backgroundColor: 'white',
    })

    viewer.addModel(pdbData, 'pdb')
    viewer.setStyle({}, { stick: {} })
    viewer.zoomTo()
    viewer.render()
  }

  return <div ref={viewerRef} style={{ width: '800px', height: '600px' }} />
}

export default dynamic(() => Promise.resolve(MoleculeViewer), {
  ssr: false, // Critical: disable SSR for WebGL
})
```

### Phase 3: API Integration (Week 3)

**Setup:**
1. Install @hey-api/openapi-ts
2. Configure to generate from backend OpenAPI spec

**Implementation:**
1. Generate type-safe API client from OpenAPI spec
2. Integrate API calls in frontend
3. Connect MoleculeViewer to backend data
4. Implement error handling and loading states

**Key Files:**
- `frontend/lib/api/client.ts` - Generated API client
- `frontend/app/visualize/page.tsx` - Integration with API

**Type-Safe API Client:**
```typescript
// Generated from OpenAPI spec
import { client } from './lib/api/client'

// Fully type-safe API calls
const convertMolecule = async (smiles: string) => {
  const response = await client.POST('/api/molecules/convert', {
    body: {
      input_format: 'smiles',
      output_format: 'pdb',
      data: smiles,
    }
  })

  // TypeScript knows exact response shape
  if (response.data) {
    return response.data // MoleculeConvertResponse type
  }
  throw new Error(response.error?.message)
}
```

### Phase 4: Docker Deployment (Week 4)

**Setup:**
1. Create Dockerfile for backend (use RDKit base image)
2. Create Dockerfile for frontend (multi-stage build)
3. Create docker-compose.yml

**Implementation:**
1. Backend Dockerfile with RDKit dependencies
2. Frontend Dockerfile with production build
3. Docker Compose with health checks
4. Environment variable configuration
5. Test full stack deployment

**Key Files:**
- `backend/Dockerfile` - Backend container
- `frontend/Dockerfile` - Frontend container
- `docker-compose.yml` - Multi-service orchestration

**docker-compose.yml:**
```yaml
version: '3.8'

services:
  backend:
    build:
      context: ./backend
      dockerfile: Dockerfile
    ports:
      - "8000:8000"
    environment:
      - ENVIRONMENT=${ENVIRONMENT:-development}
      - CORS_ORIGINS=${CORS_ORIGINS:-http://localhost:3000}
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:8000/health"]
      interval: 30s
      timeout: 10s
      retries: 3
    # Only mount volumes in development
    volumes:
      - ./backend:/app
    command: uvicorn app.main:app --host 0.0.0.0 --port 8000 --reload

  frontend:
    build:
      context: ./frontend
      dockerfile: Dockerfile
    ports:
      - "3000:3000"
    environment:
      - NEXT_PUBLIC_API_URL=${NEXT_PUBLIC_API_URL:-http://localhost:8000}
    depends_on:
      backend:
        condition: service_healthy
    healthcheck:
      test: ["CMD", "curl", "-f", "http://localhost:3000"]
      interval: 30s
      timeout: 10s
      retries: 3

# For production deployment, create docker-compose.prod.yml without volumes
```

### Phase 5: Advanced Features (Week 5)

**Caching:**
- Redis for molecule conversion results
- Browser caching for 3D models

**Performance:**
- WebGL Level of Detail (LOD) rendering
- Progressive loading for large molecules
- Web Workers for heavy calculations

**Features:**
- Multiple visualization styles
- Color schemes (by element, property)
- Measurement tools (distance, angle)
- Screenshot/export functionality

## Validation Gates

### Phase 1: Backend Foundation

```bash
# Install dependencies
cd backend
uv venv
uv sync --all-extras

# Run tests
uv run pytest tests/unit/ -v

# Start server
uv run uvicorn app.main:app --reload

# Verify OpenAPI docs
curl http://localhost:8000/docs
# Expected: Swagger UI page with API documentation

# Test SMILES validation
curl -X POST http://localhost:8000/api/molecules/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "c1ccccc1"}'
# Expected: {"valid": true, "smiles": "c1ccccc1"}

# Test invalid SMILES (fast-failure)
curl -X POST http://localhost:8000/api/molecules/validate \
  -H "Content-Type: application/json" \
  -d '{"smiles": "invalid"}'
# Expected: 422 error with actionable message

# Test format conversion
curl -X POST http://localhost:8000/api/molecules/convert \
  -H "Content-Type: application/json" \
  -d '{
    "input_format": "smiles",
    "output_format": "pdb",
    "data": "c1ccccc1"
  }'
# Expected: PDB data for benzene with molecule info
```

### Phase 2: Frontend Foundation

```bash
cd frontend
npm install
npm run dev
# Expected: Dev server starts on http://localhost:3000

# Build test
npm run build
# Expected: Production build succeeds with no TypeScript errors

# Type checking
npm run type-check
# Expected: No type errors

# Browser test (manual):
# 1. Navigate to http://localhost:3000/visualize
# 2. Should see 3D viewer with sample molecule
# 3. Test rotation, zoom, pan controls
# 4. Verify smooth 60 FPS rendering
```

### Phase 3: API Integration

```bash
# Generate API client
cd frontend
npx @hey-api/openapi-ts -i http://localhost:8000/openapi.json -o lib/api
# Expected: Type-safe client generated in lib/api/

# TypeScript compilation
npm run build
# Expected: No type errors, successful build

# End-to-end test (manual):
# 1. Start backend: cd backend && uv run uvicorn app.main:app
# 2. Start frontend: cd frontend && npm run dev
# 3. Navigate to visualizer
# 4. Enter SMILES: "c1ccccc1"
# 5. Click "Visualize"
# Expected: Benzene molecule renders in 3D within 2 seconds
```

### Phase 4: Docker Deployment

```bash
# Build and start all services
docker-compose up --build

# Check service health
docker-compose ps
# Expected: All containers running and healthy

# Test backend through Docker
curl http://localhost:8000/docs
# Expected: OpenAPI documentation accessible

# Test frontend through Docker
curl http://localhost:3000
# Expected: HTML page returned

# Test end-to-end through Docker
# 1. Navigate to http://localhost:3000
# 2. Upload SMILES or PDB file
# 3. Verify 3D visualization works
# Expected: Full functionality through containerized services

# Clean up
docker-compose down
```

### Phase 5: Advanced Features

```bash
# Performance test - Large molecule (500+ atoms)
curl -X POST http://localhost:8000/api/molecules/convert \
  -H "Content-Type: application/json" \
  -d '{
    "input_format": "smiles",
    "output_format": "pdb",
    "data": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"  # Ibuprofen
  }'
# Expected: Response < 500ms, smooth 3D rendering

# Load test - Concurrent requests
ab -n 100 -c 10 -p smiles.json -T application/json \
  http://localhost:8000/api/molecules/validate
# Expected: All requests succeed, avg response < 200ms

# Cache test (if Redis implemented)
# 1. Convert SMILES to PDB
# 2. Repeat same request
# Expected: Second request significantly faster (cached)
```

## Acceptance Criteria

### Backend Criteria

- [ ] FastAPI server running with auto-generated OpenAPI docs at `/docs`
- [ ] RDKit successfully converts between SMILES, PDB, SDF, MOL2, XYZ formats
- [ ] Pydantic v2 models validate all input/output data with type hints
- [ ] All endpoints return Chemical JSON format for visualization
- [ ] Test coverage ≥ 80% for molecule operations module
- [ ] API response time < 200ms for simple molecules (< 50 atoms)
- [ ] API response time < 500ms for complex molecules (50-500 atoms)
- [ ] Proper error handling with actionable messages (no fishy fallbacks!)
- [ ] All errors include troubleshooting guidance
- [ ] CORS configured correctly for frontend origin
- [ ] Health check endpoint returns service status

### Frontend Criteria

- [ ] Next.js app builds without TypeScript errors
- [ ] 3Dmol.js successfully renders molecules in 3D
- [ ] Interactive controls work smoothly (rotate, zoom, pan)
- [ ] Type-safe API client with full TypeScript support
- [ ] Client-side only rendering for WebGL components (no SSR errors)
- [ ] Responsive design works on desktop (1920x1080, 1366x768)
- [ ] Responsive design works on mobile (375x667, 390x844)
- [ ] Loading states display during API calls
- [ ] Error boundaries catch and display errors gracefully
- [ ] Renders at 60 FPS for molecules with 500+ atoms
- [ ] Supports visualization styles: stick, sphere, cartoon
- [ ] Color schemes work: by element, by atom type

### Integration Criteria

- [ ] Upload SMILES string → View 3D structure (< 2 seconds total)
- [ ] Upload PDB file → View 3D structure (< 3 seconds)
- [ ] Convert between formats through UI successfully
- [ ] Real-time molecule property display (formula, MW, atoms, bonds)
- [ ] Error messages from backend display in frontend
- [ ] API type mismatch caught at build time (not runtime)
- [ ] OpenAPI spec changes auto-regenerate types

### Deployment Criteria

- [ ] `docker-compose up` starts all services successfully
- [ ] Environment variables properly configured (`.env.example` provided)
- [ ] Health checks pass for all containers
- [ ] Frontend can communicate with backend through Docker network
- [ ] Logs output to stdout/stderr for monitoring
- [ ] Production build optimized (minified, tree-shaken)
- [ ] Static assets served efficiently
- [ ] Database migrations run automatically (if DB added)

## Risk Assessment & Mitigation

### Risk 1: RDKit Dependency Complexity
**Impact:** High | **Probability:** Medium

**Issue:** RDKit has complex dependencies (C++ libraries) that can be difficult to install consistently across environments.

**Mitigation:**
- Use official RDKit Docker image as base (`continuumio/miniconda3`)
- Document exact version in `pyproject.toml`: `rdkit = "^2024.03.1"`
- Provide conda environment.yml as alternative
- Include troubleshooting guide in README

**Validation:**
```bash
# Verify RDKit installation
uv run python -c "from rdkit import Chem; print(Chem.__version__)"
```

### Risk 2: WebGL Performance with Large Molecules
**Impact:** High | **Probability:** Medium

**Issue:** Molecules with 1000+ atoms may cause performance degradation in browser.

**Mitigation:**
- Implement Level of Detail (LOD) rendering
- Use Web Workers for heavy calculations
- Progressive loading for large structures
- Resolution scaling for performance vs quality tradeoff
- Display warning for molecules > 1000 atoms

**Example:**
```typescript
// Progressive loading
if (numAtoms > 1000) {
  // Load backbone first
  viewer.addModel(backboneOnly, 'pdb')
  viewer.render()

  // Then add side chains
  setTimeout(() => {
    viewer.addModel(sideChains, 'pdb')
    viewer.render()
  }, 100)
}
```

### Risk 3: CORS Issues Between Services
**Impact:** Medium | **Probability:** High

**Issue:** Browser blocks requests from Next.js (localhost:3000) to FastAPI (localhost:8000).

**Mitigation:**
- Configure FastAPI CORS middleware with environment variables:
```python
from fastapi.middleware.cors import CORSMiddleware
import os

cors_origins = os.getenv("CORS_ORIGINS", "http://localhost:3000").split(",")

app.add_middleware(
    CORSMiddleware,
    allow_origins=cors_origins,  # From environment variable
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)
```
- Use Next.js rewrites for development proxy in `next.config.js`
- Document production deployment with reverse proxy (Nginx, Traefik)

### Risk 4: Type Mismatch Between Backend/Frontend
**Impact:** Medium | **Probability:** Medium

**Issue:** Backend API changes break frontend without warning until runtime.

**Mitigation:**
- Use OpenAPI spec as single source of truth
- Automated type generation with @hey-api/openapi-ts
- CI/CD pipeline validates types are in sync
- Run type generation in pre-commit hook
- Integration tests catch mismatches

**Validation:**
```bash
# In CI/CD pipeline
npm run generate-api-client
npm run type-check
# Fails build if types don't match
```

### Risk 5: Chemistry Domain Expertise Required
**Impact:** Medium | **Probability:** Low

**Issue:** Implementing chemistry features requires domain knowledge we may not have.

**Mitigation:**
- Follow established patterns from PubChem/ChEMBL APIs
- Use RDKit Cookbook examples as reference
- Reference academic papers for visualization best practices
- Start with simple features (SMILES → 3D)
- Incremental complexity (properties → search → advanced)

**References:**
- RDKit Cookbook: https://www.rdkit.org/docs/Cookbook.html
- PubChem API: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- 3Dmol.js Examples: https://3dmol.csb.pitt.edu/doc/tutorial-code.html

## Confidence Level

**7/10** for one-pass implementation success

**Reasons for Confidence:**
- ✅ Tech stack is well-established and documented
- ✅ Clear examples exist for FastAPI + RDKit integration
- ✅ 3Dmol.js has comprehensive documentation and examples
- ✅ Next.js + WebGL patterns are proven
- ✅ Docker deployment is straightforward
- ✅ OpenAPI type generation is automated

**Risk Factors:**
- ⚠️ RDKit installation complexity
- ⚠️ WebGL performance tuning may need iteration
- ⚠️ Chemistry domain knowledge gaps

**Recommendation:**
Start with Phase 1 (Backend Foundation) to validate RDKit integration early. This de-risks the most uncertain dependency.

## References

### Official Documentation

**Frontend:**
- Next.js: https://nextjs.org/docs
- 3Dmol.js: https://3dmol.csb.pitt.edu/
- 3Dmol.js GitHub: https://github.com/3dmol/3Dmol.js
- molecule-3d-for-react: https://github.com/Autodesk/molecule-3d-for-react

**Backend:**
- FastAPI: https://fastapi.tiangolo.com/
- RDKit Documentation: https://www.rdkit.org/docs/
- RDKit Getting Started: https://www.rdkit.org/docs/GettingStartedInPython.html
- RDKit Cookbook: https://www.rdkit.org/docs/Cookbook.html
- Pydantic v2: https://docs.pydantic.dev/latest/

**API Design:**
- PubChem REST API: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest
- ChEMBL API: https://www.ebi.ac.uk/chembl/api/data/docs

### Academic Papers

- 3Dmol.js Paper: https://academic.oup.com/bioinformatics/article/31/8/1322/213186
- NGL Viewer Paper: https://academic.oup.com/nar/article/43/W1/W576/2467902
- Chemical JSON Format: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-017-0241-z

### Tutorials & Examples

- Vercel: Building WebGL with Next.js: https://vercel.com/blog/building-an-interactive-webgl-experience-in-next-js
- FastAPI + Next.js Monorepo: https://www.vintasoftware.com/blog/nextjs-fastapi-monorepo
- FastAPI Best Practices: https://github.com/zhanymkanov/fastapi-best-practices
- RDKit for AI: https://zoehlerbz.medium.com/manipulation-of-molecules-with-rdkit-in-python-for-ai-models-8023f1e677c7

### Performance & Optimization

- WebGL Best Practices: https://developer.mozilla.org/en-US/docs/Web/API/WebGL_API/WebGL_best_practices
- Next.js Performance: https://stackoverflow.blog/2022/12/20/best-practices-to-increase-the-speed-for-next-js-apps/

### Deployment

- FastAPI Docker: https://fastapi.tiangolo.com/deployment/docker/
- Next.js Docker: https://nextjs.org/docs/deployment
- Railway Deployment Guide: https://docs.railway.app/
- Render Deployment Guide: https://render.com/docs

### Chemistry Data Formats

- Chemical File Formats: https://en.wikipedia.org/wiki/Chemical_file_format
- SMILES Tutorial: https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
- PDB Format: https://www.wwpdb.org/documentation/file-format

---

## Peer Review History

### Review #1 - 2025-10-05T23:00:00Z

**Reviewer:** Context-Naive Peer Review
**Status:** Improvements Applied
**Version:** 1.0 → 1.1

**Changes Applied:**
1. ✅ Fixed PRP ID mismatch (PRP-001 → PRP-0)
2. ✅ Updated Pydantic v2 syntax (`@validator` → `@field_validator`)
3. ✅ Added validation for VisualizationPrepareRequest (at least one input required)
4. ✅ Fixed TypeScript union type syntax (Python 3.9 compatibility)
5. ✅ Added health endpoint to implementation plan
6. ✅ Fixed Docker Compose for development/production separation
7. ✅ Added environment variables for CORS and API URLs
8. ✅ Fixed MoleculeViewer component with proper Props interface and dependencies
9. ✅ Added test structure with conftest.py mention
10. ✅ Fixed UV sync command (`--all-extras` flag)

**Remaining Recommendations:**
- Consider adding example pytest fixtures in conftest.py
- Add docker-compose.prod.yml example without volumes
- Include example .env files for both services

**Quality Score:** 9/10 (excellent foundation, minor enhancements possible)

---

**Next Steps After PRP Approval:**

1. Review and approve this PRP: `/peer-review PRP-0`
2. Execute implementation: `/execute-prp PRP-0-chemistry-3d-platform-stack.md`
3. Update context after implementation: `/update-context --prp PRP-0-chemistry-3d-platform-stack.md`
4. Validate system integrity: `/validate-prp-system`
