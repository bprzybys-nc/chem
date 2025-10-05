# Suggested Commands

## Backend (Python + FastAPI)

### Setup
```bash
cd backend
uv sync --all-extras        # Install all dependencies
```

### Development
```bash
uv run uvicorn app.main:app --reload    # Start dev server
open http://localhost:8000/docs         # View API docs
```

### Testing
```bash
uv run pytest tests/ -v                 # All tests
uv run pytest tests/unit/ -v            # Unit tests only
uv run pytest tests/ --cov=app --cov-report=html  # With coverage
```

### Code Quality
```bash
uv run ruff check . --fix               # Lint and fix
uv run ruff format .                    # Format code
```

## Frontend (Next.js + TypeScript)

### Setup
```bash
cd frontend
npm install                             # Install dependencies
```

### Development
```bash
npm run dev                             # Start dev server (turbopack)
open http://localhost:3000              # View app
```

### Code Quality
```bash
npm run lint                            # ESLint check
npm run build                           # Production build test
```

### API Client Generation
```bash
npm run generate:api                    # Generate API client from openapi.json
```

## Docker (Full Stack)

```bash
docker-compose up --build               # Build and start all services
docker-compose up                       # Start services
docker-compose down                     # Stop services
```

## System Utilities (macOS Darwin)

```bash
ls -la                                  # List files with details
find . -name "*.py"                     # Find Python files
grep -r "pattern" .                     # Search for pattern
git status                              # Git status
git log --oneline -10                   # Recent commits
```