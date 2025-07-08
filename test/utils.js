// test/utils.js - Test utilities and helper functions
const fs = require('fs');
const path = require('path');
const { TEST_CONFIG, TEST_ENDPOINTS, TEST_FIXTURES } = require('./config');

function sleep(ms) {
  return new Promise(resolve => setTimeout(resolve, ms));
}

function generateTestData() {
  return {
    randomSmiles: () => {
      const smiles = TEST_FIXTURES.molecules.simple;
      return smiles[Math.floor(Math.random() * smiles.length)];
    },
    
    randomObject: () => {
      const objects = TEST_FIXTURES.objects.common;
      return objects[Math.floor(Math.random() * objects.length)];
    },
    
    randomId: () => {
      return Math.random().toString(36).substr(2, 9);
    },
    
    timestamp: () => {
      return new Date().toISOString();
    }
  };
}

class MockHttpClient {
  constructor(baseUrl = TEST_ENDPOINTS.base) {
    this.baseUrl = baseUrl;
    this.responses = new Map();
    this.requests = [];
  }
  
  mockResponse(endpoint, response, status = 200) {
    this.responses.set(endpoint, { response, status });
  }
  
  clearMocks() {
    this.responses.clear();
    this.requests = [];
  }
  
  getRequests() {
    return [...this.requests];
  }
  
  async fetch(endpoint, options = {}) {
    const url = `${this.baseUrl}${endpoint}`;
    const request = {
      url,
      endpoint,
      method: options.method || 'GET',
      body: options.body,
      headers: options.headers,
      timestamp: new Date().toISOString()
    };
    
    this.requests.push(request);
    
    if (this.responses.has(endpoint)) {
      const mock = this.responses.get(endpoint);
      await sleep(TEST_CONFIG.mock.openai.delay);
      
      return {
        ok: mock.status >= 200 && mock.status < 300,
        status: mock.status,
        json: async () => mock.response,
        text: async () => JSON.stringify(mock.response)
      };
    }
    
    throw new Error(`No mock response configured for ${endpoint}`);
  }
}

class TestFileManager {
  constructor() {
    this.tempFiles = [];
    this.testDir = TEST_CONFIG.directories.temp;
  }
  
  createTempFile(filename, content) {
    if (!fs.existsSync(this.testDir)) {
      fs.mkdirSync(this.testDir, { recursive: true });
    }
    
    const filepath = path.join(this.testDir, filename);
    fs.writeFileSync(filepath, content);
    this.tempFiles.push(filepath);
    return filepath;
  }
  
  createTestSdf(smiles, filename) {
    const sdfContent = `
  Mrv2114 12152024

  1  0  0  0  0  0            999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
M  END
> <SMILES>
${smiles}

> <_test_file>
true

$$$$
`;
    return this.createTempFile(filename || `${smiles.replace(/[^a-zA-Z0-9]/g, '_')}.sdf`, sdfContent);
  }
  
  createTestImage(imageData, filename) {
    const buffer = Buffer.from(imageData, 'base64');
    const filepath = path.join(this.testDir, filename);
    fs.writeFileSync(filepath, buffer);
    this.tempFiles.push(filepath);
    return filepath;
  }
  
  cleanup() {
    this.tempFiles.forEach(filepath => {
      if (fs.existsSync(filepath)) {
        fs.unlinkSync(filepath);
      }
    });
    this.tempFiles = [];
  }
}

class TestAssertions {
  static isValidSmiles(smiles) {
    if (typeof smiles !== 'string' || smiles.length === 0) {
      return false;
    }
    
    if (smiles.includes(' ') || smiles === 'invalid') {
      return false;
    }
    
    return /^[A-Za-z0-9\[\]()=#+\-\.@:\/\\%]+$/.test(smiles);
  }
  
  static isValidSdfPath(path) {
    return typeof path === 'string' && 
           path.endsWith('.sdf') && 
           (path.startsWith('/') || path.startsWith('http'));
  }
  
  static hasValidTestResponse(response) {
    return response && 
           response.output && 
           response.output._test &&
           response.output._test.timestamp;
  }
  
  static arrayContainsValidSmiles(array) {
    return Array.isArray(array) && 
           array.length > 0 && 
           array.every(smiles => this.isValidSmiles(smiles));
  }
  
  static responseMatchesFixture(response, fixture) {
    if (!response || !response.output) return false;
    
    const { smiles } = response.output;
    const expectedSmiles = fixture.smiles;
    
    return expectedSmiles.some(expected => 
      smiles.includes(expected)
    );
  }
}

class TestDataBuilder {
  constructor() {
    this.data = {};
  }
  
  withObject(object) {
    this.data.object = object;
    return this;
  }
  
  withSmiles(smiles) {
    this.data.smiles = Array.isArray(smiles) ? smiles : [smiles];
    return this;
  }
  
  withImage(imageBase64) {
    this.data.imageBase64 = imageBase64;
    return this;
  }
  
  withCroppedImage(croppedImageBase64) {
    this.data.croppedImageBase64 = croppedImageBase64;
    return this;
  }
  
  withCoordinates(x, y) {
    this.data.x = x;
    this.data.y = y;
    return this;
  }
  
  withTestObject(testObject) {
    this.data.testObject = testObject;
    return this;
  }
  
  build() {
    return { ...this.data };
  }
}

class TestServerUtils {
  static async waitForServer(port, timeout = 5000) {
    const start = Date.now();
    
    while (Date.now() - start < timeout) {
      try {
        const response = await fetch(`http://localhost:${port}/test/health`);
        if (response.ok) {
          return true;
        }
      } catch (error) {
        // Server not ready yet
      }
      
      await sleep(100);
    }
    
    return false;
  }
  
  static async resetTestEnvironment(port = TEST_CONFIG.server.port) {
    try {
      const response = await fetch(`http://localhost:${port}/test/utils/reset`);
      return response.ok;
    } catch (error) {
      console.error('Failed to reset test environment:', error);
      return false;
    }
  }
  
  static async getTestFixtures(port = TEST_CONFIG.server.port) {
    try {
      const response = await fetch(`http://localhost:${port}/test/fixtures`);
      return response.ok ? await response.json() : null;
    } catch (error) {
      console.error('Failed to get test fixtures:', error);
      return null;
    }
  }
}

class TestReporter {
  constructor() {
    this.results = [];
    this.startTime = Date.now();
  }
  
  logTest(testName, status, duration, details = {}) {
    const result = {
      testName,
      status,
      duration,
      details,
      timestamp: new Date().toISOString()
    };
    
    this.results.push(result);
    console.log(`[TEST] ${testName}: ${status} (${duration}ms)`);
  }
  
  generateReport() {
    const endTime = Date.now();
    const totalDuration = endTime - this.startTime;
    
    const passed = this.results.filter(r => r.status === 'PASSED').length;
    const failed = this.results.filter(r => r.status === 'FAILED').length;
    const total = this.results.length;
    
    return {
      summary: {
        total,
        passed,
        failed,
        duration: totalDuration,
        passRate: total > 0 ? ((passed / total) * 100).toFixed(1) : '0'
      },
      results: this.results
    };
  }
}

module.exports = {
  sleep,
  generateTestData,
  MockHttpClient,
  TestFileManager,
  TestAssertions,
  TestDataBuilder,
  TestServerUtils,
  TestReporter
}; 