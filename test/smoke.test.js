// Simple smoke tests - validate basic app functionality
const fs = require('fs');
const path = require('path');

describe('Basic App Smoke Tests', () => {
  test('Core files exist', () => {
    const coreFiles = [
      'frontend/core/index.html',
      'frontend/core/app.js', 
      'backend/api/server.js',
      'package.json'
    ];
    
    coreFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });
  });

  test('Package.json has required dependencies', () => {
    const packageJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'package.json')));
    
    const requiredDeps = ['express', 'cors', 'openai'];
    requiredDeps.forEach(dep => {
      expect(packageJson.dependencies[dep]).toBeDefined();
    });
  });

  test('Frontend assets exist', () => {
    const assetFiles = [
      'frontend/assets/style.css',
      'frontend/assets/camera.svg',
      'frontend/assets/account.svg'
    ];
    
    assetFiles.forEach(file => {
      expect(fs.existsSync(path.join(__dirname, '..', file))).toBe(true);
    });
  });

  test('Environment configuration is valid', () => {
    // Test that scripts directory is set up
    expect(fs.existsSync(path.join(__dirname, '..', 'scripts'))).toBe(true);
    expect(fs.existsSync(path.join(__dirname, '..', 'dev'))).toBe(true);
  });
}); 