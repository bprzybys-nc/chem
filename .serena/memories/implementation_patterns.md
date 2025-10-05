# Implementation Patterns

Key patterns extracted from current codebase implementation.

## Backend Patterns (Python + FastAPI + RDKit)

### Service Layer Pattern
**Class**: `MoleculeService` (backend/app/molecules/service.py)

**Structure**:
- All methods are `@staticmethod` (stateless)
- Grouped by functionality:
  - `from_*()` methods: Parse various formats
  - `to_*()` methods: Convert to various formats
  - `calculate_properties()`: Compute molecular properties
  - `convert()`: High-level conversion wrapper

**Benefits**:
- No instance state needed
- Easy to test (no mocking)
- Functional programming style
- Clear separation from API layer

### Fast-Failure Error Handling
**Pattern** (used in all `from_*()` methods):
```python
mol = Chem.MolFromSmiles(smiles)
if mol is None:
    raise ValueError(
        f"Invalid SMILES string: {smiles}. "
        f"Check syntax at: https://..."
    )
return mol
```

**Key aspects**:
- Check for `None` immediately after RDKit parsing
- Raise `ValueError` with specific error message
- Include troubleshooting guidance (URL or steps)
- No default fallbacks or silent failures

### Pydantic Schema Organization
**Pattern** (backend/app/molecules/schemas.py):
- 7 schema classes, all under 20 lines
- Separate Request/Response schemas per endpoint
- Nested models for complex data (`MoleculeInfo`)
- Field descriptions for OpenAPI docs
- Validation with `field_validator`

**Naming**:
- `{Operation}Request`: Input schema
- `{Operation}Response`: Output schema
- Descriptive field names (no abbreviations)

### API Router Pattern
**Structure** (backend/app/molecules/router.py):
- FastAPI router with `/api/molecules` prefix
- Exception handling at boundary (HTTPException)
- Depends on service layer for business logic
- Type-safe with Pydantic schemas

**Error handling**:
```python
try:
    result = service.method()
    return Response(...)
except ValueError as e:
    raise HTTPException(status_code=422, detail=str(e))
```

## Frontend Patterns (Next.js + React + TypeScript)

### Dynamic Import Pattern
**Component**: `MoleculeViewer` (frontend/components/MoleculeViewer.tsx)

**Pattern**:
```typescript
// Disable SSR for WebGL component
export default dynamic(() => Promise.resolve(MoleculeViewer), {
  ssr: false
})
```

**Why**:
- 3Dmol.js requires WebGL (browser-only)
- Next.js SSR would fail
- Dynamic import with `ssr: false` solves this

### Client-Side Library Loading
**Pattern** (MoleculeViewer):
```typescript
useEffect(() => {
  if (typeof window !== 'undefined' && !$3Dmol.current) {
    import('3dmol/build/3Dmol.js').then((module) => {
      $3Dmol.current = module.default
      initializeViewer()
    })
  }
}, [initializeViewer])
```

**Benefits**:
- Loads 3Dmol.js only on client
- Prevents SSR errors
- Lazy loading for performance

### Generated API Client Pattern
**Setup** (frontend/package.json):
```json
"generate:api": "openapi-ts -i openapi.json -o api -c @hey-api/client-fetch"
```

**Usage**:
- Never write manual API calls
- Generate from OpenAPI spec
- Full TypeScript types
- Located in `frontend/src/api/`

## Testing Patterns

### Real Functionality Testing
**Pattern** (backend/tests/unit/test_molecules.py):
- Test real RDKit behavior
- No mocks or stubs
- Use FastAPI `TestClient`
- Verify actual responses

**Structure**:
```python
class TestMoleculeValidation:
    def test_validate_valid_smiles(self, client, benzene_smiles):
        response = client.post("/api/molecules/validate", json={...})
        assert response.status_code == 200
        data = response.json()
        assert data["valid"] is True
```

### Test Organization
- Classes group related tests
- Fixtures in `conftest.py`
- Real chemistry data (benzene, caffeine SMILES)
- 12 tests covering all endpoints

## Configuration Patterns

### Pydantic Settings Pattern
**File**: `backend/app/config.py`
- Use `pydantic_settings.BaseSettings`
- Environment variables for runtime config
- No hardcoded values in business logic

### Environment-Based Config
**Docker Compose**:
```yaml
environment:
  - ENVIRONMENT=production
  - CORS_ORIGINS=http://localhost:3000
```

**Frontend**:
```yaml
environment:
  - NEXT_PUBLIC_API_URL=http://localhost:8000
```

## File Organization Patterns

### Module-Based Structure
**Backend**:
```
app/
├── main.py           # Entry point
├── config.py         # Settings
└── molecules/        # Feature module
    ├── router.py     # API endpoints
    ├── schemas.py    # Pydantic models
    └── service.py    # Business logic
```

**Not file-type structure** (avoided):
```
app/
├── routers/
├── schemas/
└── services/
```

### Component Structure
**Frontend**:
```
frontend/
├── app/              # App Router pages
├── components/       # Reusable components
└── src/api/          # Generated API client
```

## Lessons Learned

1. **Static service methods** work well for chemistry operations (no state needed)
2. **Fast-failure pattern** catches RDKit None returns immediately
3. **Troubleshooting in errors** reduces debugging time significantly
4. **Dynamic imports** essential for browser-only WebGL libraries
5. **Generated API clients** eliminate manual API coding and typing errors
6. **Module-based structure** scales better than file-type organization
7. **Real functionality tests** catch integration issues that mocks hide