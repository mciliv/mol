# Testing Guide

## Available Test Procedures

### **Quick Tests** (< 5 seconds)
```bash
npm run test:quick        # Smoke tests only
./scripts/test-quick      # Direct script
```
- ✅ Core files exist
- ✅ Dependencies available
- ✅ Basic functionality

### **Full Test Suite** (30+ seconds)
```bash
npm test                  # All tests (167 total)
npm run test:with-server  # With live server for integration
./scripts/test-with-server
```

### **Test Categories**
```bash
npm run test:unit         # Unit tests (frontend + backend)
npm run test:integration  # Integration tests (need server)  
npm run test:smoke        # Smoke tests (core validation)
```

### **Development Workflows**

#### **Test → LiveReload** (Most Common)
```bash
npm run procedure:test-livereload
./scripts/test-livereload
```
1. Runs full test suite
2. Starts development server with LiveReload
3. Integration tests can connect to running server

#### **Quick Development** (Skip Tests)
```bash
npm run procedure:quick-dev
./scripts/quick-dev
```
- Skips tests, immediate dev server
- For rapid iteration

#### **Debug Mode** (Maximum Visibility)
```bash
npm run procedure:debug
./scripts/debug-pipeline
```
- Maximum error visibility
- VS Code debugger attached
- All logs visible

## Test Structure

### **167 Total Tests**
- **Unit Tests**: 20+ (frontend/backend functions)
- **Integration Tests**: 140+ (full system, browser automation)
- **Smoke Tests**: 4 (basic validation)

### **Test Projects** (Jest Multi-Project)
- `unit-frontend`: DOM/browser tests (jsdom)
- `unit-backend`: Node.js server tests
- `integration`: Puppeteer + server tests  
- `smoke`: Fast validation tests

### **Test Dependencies**
- `jest` - Test runner
- `jsdom` - DOM simulation
- `puppeteer` - Browser automation
- `supertest` - HTTP testing
- `jest-fetch-mock` - Network mocking

## VS Code Integration

### **Tasks** (Cmd+Shift+P → "Tasks: Run Task")
- `test-livereload-procedure`
- `quick-dev-procedure`
- `debug-pipeline-procedure`

### **Launch Configurations** (F5)
- Debug individual components
- Full pipeline debugging

## Troubleshooting

### **Common Issues**
- **Port conflicts**: Scripts auto-cleanup ports 8080, 35729, 35730
- **Server not running**: Use `test:with-server` for integration tests
- **Dependencies missing**: `npm install` to restore

### **Error Patterns**
- `Connection refused` → Need server running
- `TextEncoder not defined` → Fixed in test setup
- `Window not defined` → Environment mismatch (fixed)

### **Quick Fixes**
```bash
# Clean everything
pkill -f "node.*server.js" && pkill -f "nodemon"

# Full dependency refresh  
rm -rf node_modules && npm install

# Test specific category
npm run test:smoke  # Just basics
``` 