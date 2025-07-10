# Test Pipeline Documentation

## Overview

The molecular visualization app uses a comprehensive 3-tier test pipeline to ensure code quality and stability throughout the development lifecycle.

## Pipeline Structure

### 1. Quick Tests (Pre-Development)
**Purpose**: Fast validation before starting development  
**Runtime**: < 5 seconds  
**Command**: `npm run test:quick`

**What it tests**:
- Server loads without crashing
- Core files are present
- Essential dependencies available
- Python/RDKit environment ready
- API endpoints respond
- Stable components unchanged

**When to use**:
- Before starting `npm run dev` (automatically runs)
- Before making changes to core components
- Quick sanity check during development

### 2. Background Tests (During Development)
**Purpose**: Comprehensive integration testing during active development  
**Runtime**: 30-60 seconds  
**Command**: `npm run test:background`

**What it tests**:
- Complete SMILES→SDF pipeline
- Schema validation coverage
- File system operations
- Error handling scenarios
- Test utilities validation
- API endpoint integration

**When to use**:
- Runs automatically in background during development
- After making significant changes
- Before committing code changes

### 3. Deployment Tests (Pre-Deployment)
**Purpose**: Full system validation before production deployment  
**Runtime**: 2-5 minutes  
**Command**: `npm run test:deployment`

**What it tests**:
- System requirements validation
- Performance under load
- Security validation
- Error recovery
- Deployment readiness
- Git repository state

**When to use**:
- Before any deployment (automatically runs)
- Before major releases
- Full system health check

## Command Reference

### Development Commands
```bash
# Safe development start (with quick tests)
npm run dev

# Unsafe development start (skip tests) 
npm run dev:unsafe

# Run background tests manually
npm run test:background
```

### Testing Commands
```bash
# Quick smoke tests only
npm run test:quick

# Standard test pipeline
npm run test:pipeline

# Full test suite including deployment
npm run test:full

# Pre-deployment validation (includes Python tests)
npm run test:pre-deploy

# Individual test categories
npm run test:unit
npm run test:integration
npm run test:fixtures
```

### Deployment Commands
```bash
# Deploy to Google Cloud (with full validation)
npm run deploy:gcp

# Quick deploy (with full validation)
npm run deploy:now

# Deploy to AWS (with full validation)
npm run deploy:aws
```

## Test File Structure

```
tests/
├── smoke.test.js      # Quick pre-development tests
├── background.test.js # Integration tests during development  
├── deployment.test.js # Comprehensive pre-deployment tests
├── server.test.js     # Core server unit tests
├── integration.test.js # Existing integration tests
├── fixtures.js        # Test data and mock objects
├── utils.js          # Test utilities and helpers
├── config.js         # Test configuration
└── setup.js          # Jest test setup
```

## Pipeline Integration

### With Development Workflow
1. **Start Development**: `npm run dev` runs quick tests first
2. **Active Development**: Background tests run automatically
3. **Code Changes**: Quick tests validate before file changes
4. **Commit Prep**: Run `npm run test:pipeline` before commits

### With Deployment Workflow
1. **Pre-Deploy**: `npm run test:pre-deploy` validates entire system
2. **Deploy**: All deploy commands include automatic validation
3. **Post-Deploy**: Can run quick smoke tests against live system

## Error Handling

### Test Failures
- **Quick tests fail**: Development blocked until fixed
- **Background tests fail**: Warning but development continues
- **Deployment tests fail**: Deployment blocked until resolved

### Recovery Commands
```bash
# If tests are stuck, clean up processes
pkill -f "jest" && npm run test:quick

# If server won't start due to port conflicts
pkill -f "node.*server.js" && lsof -ti:8080 | xargs kill -9

# Reset test environment
rm -rf node_modules/.cache/jest && npm run test:quick
```

## Performance Targets

| Test Type | Target Runtime | Max Timeout |
|-----------|---------------|-------------|
| Quick | < 5 seconds | 10 seconds |
| Background | < 60 seconds | 2 minutes |
| Deployment | < 5 minutes | 10 minutes |

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
1. **Quick tests**: Focus on basic health checks only
2. **Background tests**: Test real integration scenarios
3. **Deployment tests**: Validate production readiness

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
npm run test:quick
```

**Python/RDKit not found**:
```bash
python --version
python -c "from rdkit import Chem"
```

**Jest hanging**:
```bash
npm run test:quick -- --detectOpenHandles --forceExit
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