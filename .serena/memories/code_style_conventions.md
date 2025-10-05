# Code Style & Conventions

## Python Backend

### Style Guide
- **Line Length**: 100 characters (Ruff configured)
- **Target Version**: Python 3.11
- **Linter**: Ruff with rules: E, F, I, N, W (E501 ignored)
- **Type Hints**: Full Pydantic v2 validation (no explicit mypy yet)
- **Async**: Use async/await for FastAPI endpoints

### Naming Conventions
- **Functions/Methods**: `snake_case`
- **Classes**: `PascalCase`
- **Constants**: `UPPER_SNAKE_CASE`
- **Private**: `_leading_underscore`

### Code Quality Standards (from CLAUDE.md)
- **No Fishy Fallbacks**: Fast failure, explicit logic, no silent errors
- **No Unmarked Mocks**: Production mocks must have `# FIXME:` comments
- **Real Testing**: No hardcoded success messages, test real functionality
- **KISS Principles**: Simple solutions, clear code over clever code

### File/Function Size Guidelines
- **Files**: Target 500 lines max
- **Functions**: Target 50 lines max
- **Classes**: Target 100 lines max

### Package Management
- **NEVER edit pyproject.toml directly**
- **ALWAYS use UV commands**:
  - `uv add package-name` (production)
  - `uv add --dev package-name` (development)
  - `uv sync` (install)

## TypeScript Frontend

### Style Guide
- **ESLint**: Next.js config (eslint-config-next)
- **TypeScript**: Strict type checking
- **Components**: React 19 functional components
- **Styling**: Tailwind CSS v4

### Naming Conventions
- **Components**: `PascalCase` (e.g., `MoleculeViewer.tsx`)
- **Files**: `kebab-case` or `PascalCase` for components
- **Functions**: `camelCase`
- **Types/Interfaces**: `PascalCase`

### Project Patterns
- **App Router**: Next.js 15 app directory structure
- **API Integration**: Generated client from OpenAPI spec
- **No Manual API Calls**: Use generated `@hey-api/client-fetch`