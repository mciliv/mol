# Test Directory Structure

## Directory Organization

```
test/
├── fixtures/          # Test data files (SDF, images, etc.)
├── unit/             # Unit tests (individual components)
└── integration/      # Integration tests (full system)
```

## Naming Conventions

### JavaScript Tests
- **Unit tests**: `*.test.js` (e.g., `camera.test.js`, `unit.test.js`)
- **Integration tests**: `*.test.js` (e.g., `integration.test.js`, `smoke.test.js`)
- **Test utilities**: `*.js` (e.g., `test-runner.js`)

### Python Tests  
- **Unit tests**: `test_*.py` (e.g., `test_chem.py`, `test_dock.py`)

### Test Data
- **Fixtures**: `*.sdf`, `*.jpg`, `*.png` (e.g., `CCO.sdf`)

## Test Categories

### Unit Tests (`test/unit/`)
- Individual component testing
- Fast execution (< 5 seconds)
- Mocked dependencies
- Examples: `camera.test.js`, `unit.test.js`, `test_chem.py`

### Integration Tests (`test/integration/`)
- Full system testing
- Real API calls
- Database/file system operations
- Examples: `integration.test.js`, `smoke.test.js`, `system.test.js`

### Fixtures (`test/fixtures/`)
- Test data files
- Sample molecules, images
- Configuration files
- Examples: `CCO.sdf`, test images

## Running Tests

```bash
# Run all tests
./test

# Run specific test categories
npm test test/unit/
npm test test/integration/

# Run individual test
npm test test/unit/camera.test.js
``` 