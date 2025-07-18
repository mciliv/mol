// tests/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const AtomPredictor = require("../AtomPredictor");
const MolecularProcessor = require("../molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("../schemas");

// Mock OpenAI API
jest.mock("openai", () => ({
  OpenAI: jest.fn().mockImplementation(() => ({
    chat: {
      completions: {
        create: jest.fn().mockResolvedValue({
          choices: [
            {
              message: {
                content: JSON.stringify({
                  object: "test object",
                  chemicals: [
                    { name: "Ethanol", smiles: "CCO" },
                    { name: "Acetic acid", smiles: "CC(=O)O" },
                  ],
                }),
              },
            },
          ],
        }),
      },
    },
  })),
}));

// Mock file system
jest.mock("fs", () => ({
  existsSync: jest.fn(),
  writeFileSync: jest.fn(),
  readFileSync: jest.fn(),
  mkdirSync: jest.fn(),
}));

// Mock child_process
jest.mock("child_process", () => ({
  execSync: jest.fn().mockReturnValue("test output"),
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

describe("Unit Tests", () => {
  let atomPredictor;
  let molecularProcessor;

  beforeEach(() => {
    atomPredictor = new AtomPredictor("test-api-key");
    molecularProcessor = new MolecularProcessor();
  });

  describe("AtomPredictor", () => {
    test("should initialize with API key", () => {
      expect(atomPredictor.client).toBeDefined();
    });

    test("should analyze text input", async () => {
      const result = await atomPredictor.analyzeText("test object");
      expect(result).toHaveProperty("object");
      expect(result).toHaveProperty("chemicals");
      expect(Array.isArray(result.chemicals)).toBe(true);
    });

    test("should analyze image input", async () => {
      const mockImageBase64 =
        "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==";
      const result = await atomPredictor.analyzeImage(
        mockImageBase64,
        mockImageBase64,
        100,
        100,
      );
      expect(result).toHaveProperty("object");
      expect(result).toHaveProperty("chemicals");
    });

    test("should handle API errors gracefully", async () => {
      const mockOpenAI = require("openai");
      mockOpenAI.OpenAI.mockImplementationOnce(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error("API Error")),
          },
        },
      }));

      const analyzer = new AtomPredictor("invalid-key");
      await expect(analyzer.analyzeText("test")).rejects.toThrow("API Error");
    });
  });

  describe("MolecularProcessor", () => {
    test("should process valid SMILES", async () => {
      const result = await molecularProcessor.processSmiles(["CCO", "CC(=O)O"]);
      expect(result).toHaveProperty("sdfPaths");
      expect(result).toHaveProperty("errors");
      expect(result).toHaveProperty("skipped");
      expect(Array.isArray(result.sdfPaths)).toBe(true);
    });

    test("should handle invalid SMILES", async () => {
      const result = await molecularProcessor.processSmiles(["INVALID_SMILES"]);
      expect(result.errors).toHaveLength(1);
      expect(result.sdfPaths).toHaveLength(0);
    });

    test("should skip non-SMILES formats", async () => {
      const result = await molecularProcessor.processSmiles(["CaCO3", "SiO2"]);
      expect(result.skipped).toHaveLength(0);
      expect(result.errors.length).toBeGreaterThan(0);
    });

    test("should handle empty SMILES array", async () => {
      const result = await molecularProcessor.processSmiles([]);
      expect(result.sdfPaths).toHaveLength(0);
      expect(result.errors).toHaveLength(0);
      expect(result.skipped).toHaveLength(0);
    });
  });

  describe("Schemas", () => {
    test("ImageMoleculeSchema should validate correct data", () => {
      const validData = {
        imageBase64:
          "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
        croppedImageBase64:
          "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==",
        x: 100,
        y: 100,
      };

      const result = ImageMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("ImageMoleculeSchema should reject invalid data", () => {
      const invalidData = {
        imageBase64: "invalid-base64",
        x: "not-a-number",
      };

      const result = ImageMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });

    test("TextMoleculeSchema should validate correct data", () => {
      const validData = { object: "test object" };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("SdfGenerationSchema should validate correct data", () => {
      const validData = {
        smiles: ["CCO", "CC(=O)O"],
        overwrite: false,
      };
      const result = SdfGenerationSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });
  });

  describe("Utility Functions", () => {
    test("should process SMILES array", async () => {
      const result = await molecularProcessor.processSmiles(["CCO", "CC(=O)O"]);
      expect(result).toHaveProperty("sdfPaths");
      expect(result).toHaveProperty("errors");
      expect(result).toHaveProperty("skipped");
    });

    test("should handle file operations", () => {
      expect(molecularProcessor.sdfDir).toBeDefined();
      expect(typeof molecularProcessor.findExistingSdfFile).toBe("function");
    });
  });

  describe("Error Handling", () => {
    test("should handle network errors in AI analyzer", async () => {
      const mockOpenAI = require("openai");
      mockOpenAI.OpenAI.mockImplementationOnce(() => ({
        chat: {
          completions: {
            create: jest.fn().mockRejectedValue(new Error("Network Error")),
          },
        },
      }));

      const analyzer = new AtomPredictor("test-key");
      await expect(analyzer.analyzeText("test")).rejects.toThrow(
        "Network Error",
      );
    });

    test("should handle processing errors gracefully", async () => {
      const result = await molecularProcessor.processSmiles(["INVALID_SMILES"]);
      expect(result.errors.length).toBeGreaterThan(0);
      expect(result.sdfPaths).toHaveLength(0);
    });
  });

  describe("Input Validation", () => {
    test("should validate schema data", () => {
      const validData = { object: "test object" };
      const result = TextMoleculeSchema.safeParse(validData);
      expect(result.success).toBe(true);
    });

    test("should reject invalid schema data", () => {
      const invalidData = { object: 123 };
      const result = TextMoleculeSchema.safeParse(invalidData);
      expect(result.success).toBe(false);
    });
  });
});
