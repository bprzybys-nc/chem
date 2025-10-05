# Context Engineering Framework

This directory contains the Context Engineering (CE) framework for systematic development.

## Structure

- **PRPs/**: Product Requirements Prompts - Detailed implementation blueprints
- **examples/**: Code patterns and best practices from the project
- **templates/**: Templates for creating INITIAL feature descriptions

## Workflow

1. **Define Feature**: Create INITIAL.md with comprehensive context
2. **Generate PRP**: Use `/generate-prp` to create detailed blueprint
3. **Peer Review**: Use `/peer-review` to validate PRP quality
4. **Execute**: Use `/execute-prp` to implement feature
5. **Review Execution**: Use `/peer-review {prp} exe` to validate implementation
6. **Update Context**: Use `/update-context` to sync with codebase

## Commands

- `/generate-prp <feature>` - Create PRP from feature description
- `/execute-prp <prp-file>` - Execute PRP implementation
- `/peer-review [prp] [exe]` - Review PRP document or execution
- `/update-context` - Sync CE with codebase
- `/validate-prp-system` - Validate PRP system integrity
- `/analyze-codebase` - Analyze codebase structure and patterns

## PRP Lifecycle

```
new → executed → validated → archived
```

- **new**: PRP created, not yet implemented
- **executed**: Implementation completed, not yet validated
- **validated**: Implementation validated and working
- **archived**: Superseded or deprecated

## Best Practices

1. **Always start with context** - Provide comprehensive background
2. **Include validation gates** - Make success criteria executable
3. **Reference existing patterns** - Build on what works
4. **Keep PRPs focused** - One clear purpose per PRP
5. **Update context regularly** - Keep CE in sync with codebase
