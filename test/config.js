const path = require('path');

const TEST_CONFIG = {
  server: {
    port: process.env.TEST_PORT || 3333,
    host: 'localhost',
    timeout: 30000
  },
  
  directories: {
    sdf: path.join(__dirname, '../test_sdf_files'),
    fixtures: path.join(__dirname, 'fixtures'),
    snapshots: path.join(__dirname, 'snapshots'),
    temp: path.join(__dirname, 'temp')
  },
  
  mock: {
    openai: {
      delay: 100,
      enable: true,
      responses: {
        timeout: 5000,
        maxRetries: 3
      }
    }
  },
  
  storage: {
    type: 'memory',
    resetBetweenTests: true
  },
  
  coverage: {
    enabled: process.env.NODE_ENV === 'test',
    threshold: {
      statements: 80,
      branches: 75,
      functions: 80,
      lines: 80
    }
  },
  
  logging: {
    level: process.env.TEST_LOG_LEVEL || 'info',
    file: path.join(__dirname, '../test.log'),
    console: true
  }
};

const TEST_DB_CONFIG = {
  host: process.env.TEST_DB_HOST || 'localhost',
  port: process.env.TEST_DB_PORT || 5432,
  database: process.env.TEST_DB_NAME || 'mol_test',
  username: process.env.TEST_DB_USER || 'test_user',
  password: process.env.TEST_DB_PASS || 'test_pass'
};

const TEST_ENDPOINTS = {
  base: `http://localhost:${TEST_CONFIG.server.port}`,
  test: {
    info: '/test/info',
    fixtures: '/test/fixtures',
    imageMolecules: '/test/image-molecules',
    objectMolecules: '/test/object-molecules',
    generateSdfs: '/test/generate-sdfs',
    health: '/test/health',
    reset: '/test/utils/reset'
  },
  main: {
    imageMolecules: '/image-molecules',
    objectMolecules: '/object-molecules'
  }
};

const TEST_FIXTURES = {
  molecules: {
    simple: ['O', 'CCO', 'C'],
    complex: [
      'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
      'CC(=O)OC1=CC=CC=C1C(=O)O',
      'C1=CC=C(C=C1)C2=CC=C(C=C2)C3=CC=C(C=C3)C4=CC=C(C=C4)C5=CC=C(C=C5)C6=CC=C(C=C6)C'
    ]
  },
  objects: {
    common: ['coffee', 'wine', 'aspirin tablet', 'water', 'plastic bottle'],
    specialized: ['pharmaceutical compound', 'organic polymer', 'metal complex']
  },
  images: {
    base64: {
      blackSquare: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==',
      whiteSquare: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk+M9QDwADhgGAWjR9awAAAABJRU5ErkJggg=='
    }
  }
};

function getTestUrl(endpoint) {
  return `${TEST_ENDPOINTS.base}${endpoint}`;
}

function setupTestEnvironment() {
  const fs = require('fs');
  
  Object.values(TEST_CONFIG.directories).forEach(dir => {
    if (!fs.existsSync(dir)) {
      fs.mkdirSync(dir, { recursive: true });
    }
  });
  
  process.env.NODE_ENV = 'test';
  process.env.TEST_PORT = TEST_CONFIG.server.port;
  
  return true;
}

function cleanupTestEnvironment() {
  const fs = require('fs');
  
  const tempDir = TEST_CONFIG.directories.temp;
  if (fs.existsSync(tempDir)) {
    fs.rmSync(tempDir, { recursive: true });
    fs.mkdirSync(tempDir, { recursive: true });
  }
  
  if (TEST_CONFIG.storage.resetBetweenTests) {
    const sdfDir = TEST_CONFIG.directories.sdf;
    if (fs.existsSync(sdfDir)) {
      fs.rmSync(sdfDir, { recursive: true });
      fs.mkdirSync(sdfDir, { recursive: true });
    }
  }
  
  return true;
}

module.exports = {
  TEST_CONFIG,
  TEST_DB_CONFIG,
  TEST_ENDPOINTS,
  TEST_FIXTURES,
  getTestUrl,
  setupTestEnvironment,
  cleanupTestEnvironment
}; 