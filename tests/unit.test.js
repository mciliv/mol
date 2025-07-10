// tests/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

const request = require("supertest");
const fs = require("fs");
const path = require("path");

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

describe("Unit Tests - Component Validation", () => {

  describe("Core Server Health", () => {
    it("should load server without crashing", () => {
      expect(app).toBeDefined();
    });

    it("should serve static files", async () => {
      const response = await request(app).get("/");
      expect(response.status).toBe(200);
      expect(response.text.toLowerCase()).toContain("<!doctype html>");
    }, 3000);
  });

  describe("Essential Files Present", () => {
    it("should have required core files", () => {
      const requiredFiles = [
        "server.js",
        "index.html", 
        "app.js",
        "style.css",
        "schemas.js",
        "sdf.py"
      ];
      
      requiredFiles.forEach(file => {
        expect(fs.existsSync(path.join(__dirname, "..", file))).toBe(true);
      });
    });

    it("should have stable components unchanged", () => {
      // Validate core SMILES engine integrity
      const sdfPy = fs.readFileSync(path.join(__dirname, "..", "sdf.py"), "utf8");
      expect(sdfPy).toContain("from rdkit import Chem");
      expect(sdfPy).toContain("def sdf");
      
      const schemas = fs.readFileSync(path.join(__dirname, "..", "schemas.js"), "utf8");
      expect(schemas).toContain("SdfGenerationSchema");
    });
  });

  describe("API Endpoints Responding", () => {
    it("should respond to generate-sdfs endpoint", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: [] });
      
      // Should respond (even with error for empty array)
      expect([200, 400, 500]).toContain(response.status);
    }, 2000);

    it("should respond to object-molecules endpoint", async () => {
      const response = await request(app)
        .post("/object-molecules")
        .send({ object: "test" });
      
      // Should respond (may fail due to missing API key, but server should handle)
      expect([200, 400, 500]).toContain(response.status);
    }, 2000);
  });

  describe("Python Dependencies", () => {
    it("should have Python available", () => {
      const { execSync } = require("child_process");
      try {
        const pythonVersion = execSync("python --version", { encoding: "utf8" });
        expect(pythonVersion).toContain("Python");
      } catch (error) {
        throw new Error("Python not available in PATH");
      }
    });

    it("should have RDKit available", () => {
      const { execSync } = require("child_process");
      try {
        execSync("python -c 'from rdkit import Chem; print(\"RDKit OK\")'", { encoding: "utf8" });
      } catch (error) {
        throw new Error("RDKit not available - required for SMILES processing");
      }
    });
  });
});

// Additional unit tests for API validation  
describe("Unit Tests - API Validation", () => {
  describe("POST /generate-sdfs endpoint behavior", () => {
    it("should return 400 if smiles array is missing", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ invalid: "data" });

      expect(response.status).toBe(400);
      expect(response.body.error).toBe("smiles array is required");
    });

    it("should handle empty smiles array", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: [] });

      expect(response.status).toBe(200);
      expect(response.body.sdfPaths).toEqual([]);
      expect(response.body.message).toContain("Generated 0 3D structures from 0 SMILES");
    });

    it("should process valid SMILES", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: ["O"] }); // Simple water molecule

      expect(response.status).toBe(200);
      expect(response.body.sdfPaths).toBeDefined();
      expect(Array.isArray(response.body.sdfPaths)).toBe(true);
      expect(response.body.message).toContain("3D structures");
    }, 10000);
  });

  describe("Schema validation components", () => {
    it("should have valid schemas defined", () => {
      const schemas = require("../schemas");
      expect(schemas.SdfGenerationSchema).toBeDefined();
      expect(schemas.ImageMoleculeSchema).toBeDefined();
      expect(schemas.TextMoleculeSchema).toBeDefined();
    });

    it("should validate SMILES array structure", () => {
      const { SdfGenerationSchema } = require("../schemas");
      
      // Valid data should parse
      expect(() => SdfGenerationSchema.parse({ smiles: ["O", "CCO"] })).not.toThrow();
      
      // Invalid data should throw
      expect(() => SdfGenerationSchema.parse({ smiles: "not an array" })).toThrow();
      expect(() => SdfGenerationSchema.parse({})).toThrow();
    });
  });
}); 