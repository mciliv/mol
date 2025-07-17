// tests/integration.test.js - Integration tests for component interactions  
// These tests validate how different parts of the system work together (30-60 seconds)

// Mock external dependencies BEFORE importing server
jest.mock('openai', () => ({
  OpenAI: jest.fn().mockImplementation(() => ({
    chat: {
      completions: {
        create: jest.fn().mockResolvedValue({
          choices: [{
            message: {
              content: JSON.stringify({
                object: "test object",
                smiles: ["CCO", "CC(=O)O"]
              })
            }
          }]
        })
      }
    }
  }))
}));

jest.mock('child_process', () => ({
  execSync: jest.fn().mockReturnValue('test output'),
  spawn: jest.fn().mockImplementation((command, args) => {
    const EventEmitter = require('events');
    const mockProcess = new EventEmitter();
    const fs = require('fs');
    const path = require('path');
    
    // Extract SMILES from args (assuming it's the second argument)
    const smiles = args[0];
    const sdfDir = args[2]; // --dir argument
    
    // Create a mock SDF file
    const sdfContent = `${smiles}
  RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$`;
    
    const sdfPath = path.join(sdfDir, `${smiles}.sdf`);
    fs.writeFileSync(sdfPath, sdfContent);
    
    // Simulate successful SDF generation
    setTimeout(() => {
      mockProcess.emit('close', 0);
    }, 10);
    
    return mockProcess;
  })
}));

const request = require('supertest');
const app = require('../server');
const fs = require('fs');
const path = require('path');

describe('Integration Tests', () => {
  beforeEach(() => {
    // Clear any test files
    const sdfDir = path.join(__dirname, '../sdf_files');
    if (fs.existsSync(sdfDir)) {
      const files = fs.readdirSync(sdfDir);
      files.forEach(file => {
        if (file.endsWith('.sdf')) {
          fs.unlinkSync(path.join(sdfDir, file));
        }
      });
    }
  });

  describe('API Endpoints', () => {
    test('POST /image-molecules should process image analysis', async () => {
      const mockImageBase64 = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==';
      
      const response = await request(app)
        .post('/image-molecules')
        .send({
          imageBase64: mockImageBase64,
          croppedImageBase64: mockImageBase64,
          x: 100,
          y: 100
        })
        .expect(200);

      expect(response.body).toHaveProperty('output');
      expect(response.body.output).toHaveProperty('object');
      expect(response.body.output).toHaveProperty('chemicals');
      expect(Array.isArray(response.body.output.chemicals)).toBe(true);
      
      // Extract SMILES from chemicals for testing
      const smiles = response.body.output.chemicals.map(chem => chem.smiles).filter(Boolean);
      expect(Array.isArray(smiles)).toBe(true);
    });

    test('POST /object-molecules should process text analysis', async () => {
      const response = await request(app)
        .post('/object-molecules')
        .send({ object: 'water bottle' })
        .expect(200);

      expect(response.body).toHaveProperty('output');
      expect(response.body.output).toHaveProperty('object');
      expect(response.body.output).toHaveProperty('chemicals');
      expect(Array.isArray(response.body.output.chemicals)).toBe(true);
      
      // Extract SMILES from chemicals for testing
      const smiles = response.body.output.chemicals.map(chem => chem.smiles).filter(Boolean);
      expect(Array.isArray(smiles)).toBe(true);
    });

    test('POST /generate-sdfs should process SMILES and generate SDF files', async () => {
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: ['CCO', 'CC(=O)O'],
          overwrite: false 
        })
        .expect(200);

      expect(response.body).toHaveProperty('sdfPaths');
      expect(response.body).toHaveProperty('errors');
      expect(response.body).toHaveProperty('skipped');
      expect(response.body).toHaveProperty('message');
      expect(Array.isArray(response.body.sdfPaths)).toBe(true);
    });

    test('GET / should serve the main HTML page', async () => {
      const response = await request(app)
        .get('/')
        .expect(200);

      expect(response.text).toContain('<!doctype html>');
      expect(response.text).toContain('Atomizer - Molecular Analysis');
      expect(response.text).toContain('app.js');
    });

    test('GET /sdf_files should serve SDF files', async () => {
      // Create a test SDF file
      const sdfDir = path.join(__dirname, '../sdf_files');
      if (!fs.existsSync(sdfDir)) {
        fs.mkdirSync(sdfDir, { recursive: true });
      }
      
      const testSdfContent = `CCO
  RDKit          3D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0  0  0  0
  2  3  1  0  0  0  0
M  END
$$$$`;
      
      fs.writeFileSync(path.join(sdfDir, 'test.sdf'), testSdfContent);

      const response = await request(app)
        .get('/sdf_files/test.sdf')
        .expect(200);

      expect(response.text).toBeDefined();
      expect(response.text).toContain('CCO');
      expect(response.text).toContain('$$$$');
    });
  });

  describe('Error Handling', () => {
    test('POST /image-molecules should handle missing image data', async () => {
      const response = await request(app)
        .post('/image-molecules')
        .send({ x: 100, y: 100 })
        .expect(400);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toBe('Invalid input data');
    });

    test('POST /object-molecules should handle missing object description', async () => {
      const response = await request(app)
        .post('/object-molecules')
        .send({})
        .expect(400);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toBe('Invalid input data');
    });

    test('POST /generate-sdfs should handle missing SMILES array', async () => {
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ overwrite: false })
        .expect(400);

      expect(response.body).toHaveProperty('error');
      expect(response.body.error).toBe('Invalid input data');
    });

    test('POST /generate-sdfs should handle invalid SMILES', async () => {
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: ['INVALID_SMILES'],
          overwrite: false 
        })
        .expect(200);

      expect(response.body.errors).toHaveLength(1);
      expect(response.body.sdfPaths).toHaveLength(0);
    });
  });

  describe('Data Flow', () => {
    test('should handle complete flow: image analysis -> SDF generation', async () => {
      const mockImageBase64 = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==';
      
      // Step 1: Analyze image
      const analysisResponse = await request(app)
        .post('/image-molecules')
        .send({
          imageBase64: mockImageBase64,
          croppedImageBase64: mockImageBase64,
          x: 100,
          y: 100
        })
        .expect(200);

      const smiles = analysisResponse.body.output.chemicals.map(chem => chem.smiles).filter(Boolean);
      expect(Array.isArray(smiles)).toBe(true);
      expect(smiles.length).toBeGreaterThan(0);

      // Step 2: Generate SDF files
      const sdfResponse = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: smiles,
          overwrite: false 
        })
        .expect(200);

      expect(sdfResponse.body.sdfPaths).toHaveLength(smiles.length);
      expect(sdfResponse.body.errors).toHaveLength(0);
    });

    test('should handle complete flow: text analysis -> SDF generation', async () => {
      // Step 1: Analyze text
      const analysisResponse = await request(app)
        .post('/object-molecules')
        .send({ object: 'water bottle' })
        .expect(200);

      const smiles = analysisResponse.body.output.chemicals.map(chem => chem.smiles).filter(Boolean);
      expect(Array.isArray(smiles)).toBe(true);

      // Step 2: Generate SDF files
      const sdfResponse = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: smiles,
          overwrite: false 
        })
        .expect(200);

      expect(Array.isArray(sdfResponse.body.sdfPaths)).toBe(true);
    });
  });

  describe('File Management', () => {
    test('should create SDF files in correct directory', async () => {
      const sdfDir = path.join(__dirname, '../sdf_files');
      
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: ['CCO'],
          overwrite: false 
        })
        .expect(200);

      expect(response.body.sdfPaths).toHaveLength(1);
      
      const sdfPath = response.body.sdfPaths[0];
      const fullPath = path.join(__dirname, '..', sdfPath);
      
      expect(fs.existsSync(fullPath)).toBe(true);
      expect(fs.readFileSync(fullPath, 'utf8')).toContain('CCO');
    });

    test('should handle overwrite flag correctly', async () => {
      const sdfDir = path.join(__dirname, '../sdf_files');
      
      // First generation
      const response1 = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: ['CCO'],
          overwrite: false 
        })
        .expect(200);

      // Second generation with overwrite
      const response2 = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: ['CCO'],
          overwrite: true 
        })
        .expect(200);

      expect(response1.body.sdfPaths).toHaveLength(1);
      expect(response2.body.sdfPaths).toHaveLength(1);
    });
  });

  describe('Schema Validation', () => {
    test('should validate image analysis request schema', async () => {
      const invalidData = {
        imageBase64: 'invalid-base64',
        x: 'not-a-number'
      };

      const response = await request(app)
        .post('/image-molecules')
        .send(invalidData)
        .expect(400);

      expect(response.body).toHaveProperty('error');
    });

    test('should validate SDF generation request schema', async () => {
      const invalidData = {
        smiles: 'not-an-array',
        overwrite: 'not-a-boolean'
      };

      const response = await request(app)
        .post('/generate-sdfs')
        .send(invalidData)
        .expect(400);

      expect(response.body).toHaveProperty('error');
    });
  });

  describe('Performance', () => {
    test('should handle multiple concurrent requests', async () => {
      const requests = Array(5).fill().map(() => 
        request(app)
          .post('/object-molecules')
          .send({ object: 'test object' })
      );

      const responses = await Promise.all(requests);
      
      responses.forEach(response => {
        expect(response.status).toBe(200);
        expect(response.body).toHaveProperty('output');
      });
    });

    test('should handle large SMILES arrays', async () => {
      const largeSmilesArray = Array(10).fill('CCO');
      
      const response = await request(app)
        .post('/generate-sdfs')
        .send({ 
          smiles: largeSmilesArray,
          overwrite: false 
        })
        .expect(200);

      expect(response.body.sdfPaths).toHaveLength(10);
    });
  });
}); 