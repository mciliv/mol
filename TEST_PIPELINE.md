# Test Pipeline Documentation

## Overview

The molecular visualization app uses a standard 3-tier test pyramid to ensure code quality and reliability at every level.

## Test Categories

### 1. Unit Tests
**Purpose**: Test individual components and functions in isolation  
**Runtime**: < 5 seconds  
**Command**: `npm run test:unit`

**What it tests**:
- Server loads without crashing
- Core files are present
- Essential dependencies available
- Python/RDKit environment ready
- API endpoints respond
- Stable components unchanged
- Individual module functionality

**When to use**:
- Before starting `npm run dev` (automatically runs)
- During active development
- Quick feedback loop
- Continuous development validation

### 2. Integration Tests  
**Purpose**: Test how different components work together
**Runtime**: 30-60 seconds  
**Command**: `npm run test:integration`

**What it tests**:
- Complete SMILES→SDF pipeline
- Schema validation coverage
- File system operations
- Error handling scenarios
- Test utilities validation
- API endpoint integration
- Multi-component workflows

**When to use**:
- After significant changes
- Before committing code changes
- Validating component interactions
- Feature completion testing

### 3. System Tests
**Purpose**: End-to-end testing of the complete system
**Runtime**: 2-5 minutes  
**Command**: `npm run test:system`

**What it tests**:
- System requirements validation
- Performance under load
- Security validation
- Error recovery
- Deployment readiness
- Git repository state
- Full user workflows

**When to use**:
- Before deployment
- Major releases
- Complete system validation
- Production readiness checks

## Command Reference

### Development Commands
```bash
# Safe development start (with unit tests)
npm run dev

# Unsafe development start (skip tests) 
npm run dev:unsafe
```

### Testing Commands
```bash
# Individual test categories
npm run test:unit         # Fast unit tests (< 5s)
npm run test:integration  # Integration tests (30-60s)  
npm run test:system      # System tests (2-5min)

# Combined test runs
npm run test             # Unit + Integration (standard)
npm run test:all         # All three test levels
npm run test:pre-deploy  # All tests + Python tests

# Development utilities
npm run test:watch       # Watch mode for active development
npm run test:legacy      # Legacy server tests
npm run test:fixtures    # Test fixture validation
```

### Deployment Commands
```bash
# Deploy to Google Cloud Functions (with full validation)
npm run deploy
```

## Test File Structure

```
tests/
├── unit.test.js           # Unit tests - individual components  
├── integration.test.js    # Integration tests - component interactions
├── system.test.js         # System tests - end-to-end workflows
├── server.test.js         # Legacy server tests (test:legacy)
├── integration.existing.js # Backup of original integration tests
├── fixtures.js            # Test data and mock objects
├── utils.js               # Test utilities and helpers
├── config.js              # Test configuration
└── setup.js               # Jest test setup
```

## Pipeline Integration

### With Development Workflow
1. **Start Development**: `npm run dev` runs unit tests first
2. **Active Development**: `npm run test:watch` for continuous feedback
3. **Feature Complete**: `npm run test:integration` validates interactions  
4. **Commit Prep**: `npm run test` (unit + integration) before commits

### With Deployment Workflow
1. **Pre-Deploy**: `npm run test:pre-deploy` validates entire system
2. **Deploy**: All deploy commands include automatic validation
3. **Post-Deploy**: `npm run test:system` can validate live system

## Error Handling

### Test Failures
- **Unit tests fail**: Development blocked until fixed (fast feedback)
- **Integration tests fail**: Feature incomplete, fix component interactions
- **System tests fail**: Deployment blocked until resolved

### Recovery Commands
```bash
# If tests are stuck, clean up processes
pkill -f "jest" && npm run test:unit

# If server won't start due to port conflicts
pkill -f "node.*server.js" && lsof -ti:8080 | xargs kill -9

# Reset test environment
rm -rf node_modules/.cache/jest && npm run test:unit
```

## Performance Targets

| Test Type | Target Runtime | Max Timeout | Purpose |
|-----------|---------------|-------------|---------|
| Unit | < 5 seconds | 10 seconds | Fast feedback loop |
| Integration | < 60 seconds | 2 minutes | Component interactions |
| System | < 5 minutes | 10 minutes | End-to-end validation |

## Configuration

### Jest Configuration (package.json)
```json
{
  "jest": {
    "testEnvironment": "node",
    "setupFilesAfterEnv": ["<rootDir>/tests/setup.js"]
  }
}
```

### Test Environment Variables
- Tests run with minimal environment setup
- OpenAI API key not required for most tests
- Python environment must be available

## Best Practices

### Writing Tests
1. **Unit tests**: Test individual functions and components in isolation
2. **Integration tests**: Test real component interactions and workflows
3. **System tests**: Test complete user scenarios and production readiness

### Test Data
- Use fixtures from `tests/fixtures.js`
- Mock external APIs when possible
- Clean up test files in `afterAll()`

### Error Messages
- Provide clear failure descriptions
- Include suggested fixes where possible
- Log relevant context for debugging

## Troubleshooting

### Common Issues

**Port conflicts during tests**:
```bash
lsof -ti:8080 | xargs kill -9
npm run test:unit
```

**Python/RDKit not found**:
```bash
python --version
python -c "from rdkit import Chem"
```

**Jest hanging**:
```bash
npm run test:unit -- --detectOpenHandles --forceExit
```

**File permission errors**:
```bash
chmod +x scripts/*.sh
npm run test:quick
```

### Test Coverage

The pipeline ensures coverage of:
- ✅ All HTTP endpoints
- ✅ SMILES processing pipeline  
- ✅ Schema validation
- ✅ Error handling
- ✅ File operations
- ✅ Python integration
- ✅ Security validation
- ✅ Performance requirements

## Integration with CI/CD

This test pipeline is designed to integrate with:
- GitHub Actions
- Google Cloud Build
- AWS CodePipeline
- Local development workflows

### Example GitHub Action
```yaml
name: Test Pipeline
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-node@v3
      - run: npm install
      - run: npm run test:pre-deploy
``` 