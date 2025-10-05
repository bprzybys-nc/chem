# Task Completion Checklist

## When Backend Task is Complete

### 1. Code Quality
- [ ] Run `uv run ruff check . --fix` - Lint and auto-fix
- [ ] Run `uv run ruff format .` - Format code
- [ ] Verify no Ruff errors remain

### 2. Testing
- [ ] Run `uv run pytest tests/unit/ -v` - Unit tests pass
- [ ] Add new tests for new functionality
- [ ] Verify 12+ tests passing (or more if added)

### 3. Type Safety
- [ ] Pydantic models validate correctly
- [ ] No runtime validation errors

### 4. Fast Failure Compliance
- [ ] No fishy fallbacks or silent error masking
- [ ] All exceptions bubble up properly
- [ ] Production mocks marked with `# FIXME:`
- [ ] Error messages include troubleshooting guidance

## When Frontend Task is Complete

### 1. Code Quality
- [ ] Run `npm run lint` - ESLint check passes
- [ ] Run `npm run build` - Production build succeeds
- [ ] No TypeScript errors

### 2. API Client
- [ ] If OpenAPI changed: `npm run generate:api`
- [ ] Verify generated client types are correct

### 3. Testing
- [ ] Manual test in browser at http://localhost:3000
- [ ] Verify 3D visualization renders correctly
- [ ] Test API integration

## Full Stack Validation

### Docker Integration
- [ ] Run `docker-compose up --build`
- [ ] Backend health check passes (http://localhost:8000/health)
- [ ] Frontend loads (http://localhost:3000)
- [ ] Frontend can communicate with backend

### Documentation
- [ ] Update README.md if new features added
- [ ] Update API documentation if endpoints changed
- [ ] Commit messages follow conventions

## TDD Workflow (from CLAUDE.md)
1. Write test first - Define expected behavior
2. Watch it fail - Ensure test validates logic
3. Write minimal code - Just enough to pass
4. Refactor - Improve while keeping tests green
5. Tests must invoke actual functionality - no mocks unless explicit