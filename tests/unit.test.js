// tests/unit.test.js - Unit tests for individual components
// These tests run quickly (< 5 seconds) and validate individual functions and modules

const request = require("supertest");
const fs = require("fs");
const path = require("path");

describe("Unit Tests - Component Validation", () => {
  let app;
  
  beforeAll(() => {
    // Minimal app setup for testing
    try {
      app = require("../server");
    } catch (error) {
      console.error("Failed to load server:", error.message);
      throw error;
    }
  });

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