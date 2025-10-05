# Chemistry Backend API

FastAPI backend for 3D molecule visualization platform with RDKit integration.

## Features

- **Molecule Format Conversion**: Convert between SMILES, PDB, SDF, MOL2, XYZ formats
- **Molecule Validation**: Validate SMILES strings with canonical form
- **Property Calculation**: Calculate molecular properties (MW, formula, LogP, etc.)
- **RESTful API**: Auto-generated OpenAPI documentation
- **Type Safety**: Full Pydantic v2 validation

## Quick Start

```bash
# Install dependencies
uv sync --all-extras

# Run tests
uv run pytest tests/unit/ -v

# Start development server
uv run uvicorn app.main:app --reload

# View API documentation
open http://localhost:8000/docs
```

## API Endpoints

- `GET /health` - Health check
- `POST /api/molecules/validate` - Validate SMILES
- `POST /api/molecules/convert` - Convert molecule formats
- `POST /api/molecules/properties` - Calculate properties

## Project Structure

```
backend/
├── app/
│   ├── main.py              # FastAPI entry point
│   ├── config.py            # Configuration
│   └── molecules/           # Molecule operations module
│       ├── router.py        # API endpoints
│       ├── schemas.py       # Pydantic models
│       └── service.py       # RDKit operations
├── tests/
│   ├── conftest.py          # Pytest fixtures
│   └── unit/                # Unit tests
└── pyproject.toml           # Dependencies (UV managed)
```

## Development

### Running Tests

```bash
# All tests
uv run pytest tests/ -v

# Unit tests only
uv run pytest tests/unit/ -v

# With coverage
uv run pytest tests/ --cov=app --cov-report=html
```

### Code Quality

```bash
# Format code
uv run ruff check . --fix

# Type checking (via Pydantic models)
uv run python -m app.main
```

## Dependencies

- **FastAPI**: Async web framework
- **Pydantic v2**: Data validation
- **RDKit**: Chemistry processing
- **UV**: Package management

## Phase 1 Status

✅ FastAPI server with CORS
✅ Health check endpoint
✅ Molecule validation endpoint
✅ Format conversion (SMILES ↔ PDB)
✅ Property calculation
✅ Unit tests (12/12 passing)
✅ OpenAPI documentation
✅ Fast-failure error handling

Next: Phase 2 - Frontend Foundation (Next.js + 3Dmol.js)
