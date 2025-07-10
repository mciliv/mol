// tests/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const AIAnalyzer = require('../ai-analyzer');
const MolecularProcessor = require('../molecular-processor');
const { ImageMoleculeSchema, TextMoleculeSchema, SdfGenerationSchema } = require('../schemas');

// Mock OpenAI API
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

// Mock file system
jest.mock('fs', () => ({
  existsSync: jest.fn(),
  writeFileSync: jest.fn(),
  readFileSync: jest.fn(),
  mkdirSync: jest.fn()
}));

// Mock child_process
jest.mock('child_process', () => ({
  execSync: jest.fn().mockReturnValue('test output')
}));

// Global app instance for all tests
let app;

beforeAll(() => {
  // Load server for all tests
  try {
    app = require("../server");
  } catch (error) {
    console.error("Failed to load server:", error.message);
    throw error;
  }
});

describe('Unit Tests', () => {
  let aiAnalyzer;
  let molecularProcessor;

  beforeEach(() => {
    aiAnalyzer = new AIAnalyzer('test-api-key');
    molecularProcessor = new MolecularProcessor();
  });

  describe('AIAnalyzer', () => {
    test('should initialize with API key', () => {
      expect(aiAnalyzer.apiKey).toBe('test-api-key');
    });

    test('should analyze text input', async () => {
      const result = await aiAnalyzer.analyzeText('test object');
      expect(result).toHaveProperty('object');
      expect(result).toHaveProperty('smiles');
      expect(Array.isArray(result.smiles)).toBe(true);
    });

    test('should analyze image input', async () => {
      const mockImageBase64 = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==';
      const result = await aiAnalyzer.analyzeImage(mockImageBase64, mockImageBase64, 100, 100);
      expect(result).toHaveProperty('object');
      expect(result).toHaveProperty('smiles');
    });

    test('should handle API errors gracefully', async () => {
      const mockOpenAI = require('openai');
      mockOpenAI.OpenAI.mockImplementationOnce(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error('API Error'))
          }
        }
      }));

      const analyzer = new AIAnalyzer('invalid-key');
      await expect(analyzer.analyzeText('test')).rejects.toThrow('API Error');
    });
  });

  describe('MolecularProcessor', () => {
    test('should process valid SMILES', async () => {
      const result = await molecularProcessor.processSmiles(['CCO', 'CC(=O)O']);
      expect(result).toHaveProperty('sdfPaths');
      expect(result).toHaveProperty('errors');
      expect(result).toHaveProperty('skipped');
      expect(Array.isArray(result.sdfPaths)).toBe(true);
    });

    test('should handle invalid SMILES', async () => {
      const result = await molecularProcessor.processSmiles(['INVALID_SMILES']);
      expect(result.errors).toHaveLength(1);
      expect(result.sdfPaths).toHaveLength(0);
    });

    test('should skip non-SMILES formats', async () => {
      const result = await molecularProcessor.processSmiles(['CaCO3', 'SiO2']);
      expect(result.skipped).toHaveLength(2);
      expect(result.sdfPaths).toHaveLength(0);
    });

    test('should handle empty SMILES array', async () => {
      const result = await molecularProcessor.processSmiles([]);
      expect(result.sdfPaths).toHaveLength(0);
      expect(result.errors).toHaveLength(0);
      expect(result.skipped).toHaveLength(0);
    });
  });

  describe('Schemas', () => {
    test('ImageMoleculeSchema should validate correct data', () => {
      const validData = {
        imageBase64: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==',
        croppedImageBase64: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==',
        x: 100,
        y: 100
      };
      
      const result = ImageMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test('ImageMoleculeSchema should reject invalid data', () => {
      const invalidData = {
        imageBase64: 'invalid-base64',
        x: 'not-a-number'
      };
      
      const result = ImageMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });

    test('TextMoleculeSchema should validate correct data', () => {
      const validData = { object: 'test object' };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test('SdfGenerationSchema should validate correct data', () => {
      const validData = {
        smiles: ['CCO', 'CC(=O)O'],
        overwrite: false
      };
      const result = SdfGenerationSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });
  });

  describe('Utility Functions', () => {
    test('should validate SMILES format', () => {
      const validSmiles = ['CCO', 'CC(=O)O', 'C1=CC=CC=C1'];
      const invalidSmiles = ['CaCO3', 'SiO2', 'INVALID'];
      
      validSmiles.forEach(smiles => {
        expect(molecularProcessor.isValidSmiles(smiles)).toBe(true);
      });
      
      invalidSmiles.forEach(smiles => {
        expect(molecularProcessor.isValidSmiles(smiles)).toBe(false);
      });
    });

    test('should generate safe filenames', () => {
      const testCases = [
        { input: 'CCO', expected: 'CCO.sdf' },
        { input: 'CC(=O)O', expected: 'CC(=O)O.sdf' },
        { input: 'C1=CC=CC=C1', expected: 'C1=CC=CC=C1.sdf' }
      ];
      
      testCases.forEach(({ input, expected }) => {
        const filename = molecularProcessor.generateSafeFilename(input);
        expect(filename).toBe(expected);
      });
    });
  });

  describe('Error Handling', () => {
    test('should handle network errors in AI analyzer', async () => {
      const mockOpenAI = require('openai');
      mockOpenAI.OpenAI.mockImplementationOnce(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error('Network Error'))
          }
        }
      }));

      const analyzer = new AIAnalyzer('test-key');
      await expect(analyzer.analyzeText('test')).rejects.toThrow('Network Error');
    });

    test('should handle file system errors in molecular processor', async () => {
      const fs = require('fs');
      fs.writeFileSync.mockImplementationOnce(() => {
        throw new Error('File system error');
      });

      await expect(molecularProcessor.processSmiles(['CCO'])).rejects.toThrow('File system error');
    });

    test('should handle Python subprocess errors', async () => {
      const { execSync } = require('child_process');
      execSync.mockImplementationOnce(() => {
        throw new Error('Python script failed');
      });

      await expect(molecularProcessor.processSmiles(['CCO'])).rejects.toThrow('Python script failed');
    });
  });

  describe('Input Validation', () => {
    test('should validate base64 image data', () => {
      const validBase64 = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==';
      const invalidBase64 = 'not-base64-data';
      
      expect(aiAnalyzer.isValidBase64(validBase64)).toBe(true);
      expect(aiAnalyzer.isValidBase64(invalidBase64)).toBe(false);
    });

    test('should validate object descriptions', () => {
      const validObjects = ['water bottle', 'coffee cup', 'plant'];
      const invalidObjects = ['', '   ', null, undefined];
      
      validObjects.forEach(obj => {
        expect(aiAnalyzer.isValidObjectDescription(obj)).toBe(true);
      });
      
      invalidObjects.forEach(obj => {
        expect(aiAnalyzer.isValidObjectDescription(obj)).toBe(false);
      });
    });
  });
}); 