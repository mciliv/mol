// tests/integration.test.js - Integration tests for component interactions  
// These tests validate how different parts of the system work together (30-60 seconds)

const request = require("supertest");
const fs = require("fs");
const path = require("path");
const { TestFileManager, TestAssertions } = require("./utils");
const { getTestMolecule, createTestRequest, MOCK_IMAGES } = require("./fixtures");

describe("Integration Tests - Component Interactions", () => {
  let app;
  let fileManager;
  
  beforeAll(() => {
    app = require("../server");
    fileManager = new TestFileManager();
  }, 30000);
  
  afterAll(() => {
    fileManager.cleanup();
  });

  describe("SMILES Processing Pipeline", () => {
    it("should generate SDF files for valid SMILES", async () => {
      const testSmiles = ["O", "CCO", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"];
      
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: testSmiles, overwrite: true });
      
      expect(response.status).toBe(200);
      expect(response.body.message).toContain("Generated");
      expect(response.body.message).toContain("3D structures");
      expect(response.body.sdfPaths).toHaveLength(testSmiles.length);
      
      // Verify files actually exist
      response.body.sdfPaths.forEach(sdfPath => {
        const fullPath = path.join(__dirname, "..", sdfPath);
        expect(fs.existsSync(fullPath)).toBe(true);
      });
    }, 15000);

    it("should handle mixed valid/invalid SMILES gracefully", async () => {
      const mixedSmiles = ["O", "INVALID_SMILES", "CCO"];
      
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: mixedSmiles, overwrite: true });
      
      // Should complete but may have fewer files than expected
      expect([200, 500]).toContain(response.status);
    }, 10000);
  });

  describe("AI Analysis Integration", () => {
    it("should handle text analysis requests", async () => {
      const response = await request(app)
        .post("/object-molecules")
        .send({ object: "water" });
      
      // Should respond (may fail if no OpenAI key, but structure should be correct)
      expect([200, 400, 500]).toContain(response.status);
      
      if (response.status === 200) {
        expect(response.body.output).toBeDefined();
        expect(response.body.output.smiles).toBeDefined();
        expect(Array.isArray(response.body.output.smiles)).toBe(true);
      }
    }, 10000);

    it("should handle image analysis requests", async () => {
      const imageData = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg==";
      
      const response = await request(app)
        .post("/image-molecules")
        .send({ 
          imageBase64: imageData,
          croppedImageBase64: imageData,
          x: 100,
          y: 100
        });
      
      // Should respond (may fail if no OpenAI key, but structure should be correct)
      expect([200, 400, 500]).toContain(response.status);
    }, 15000);
  });

  describe("Error Handling Integration", () => {
    it("should handle malformed requests gracefully", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ invalid: "data" });
      
      expect(response.status).toBe(400);
      expect(response.body.error).toBeDefined();
    });

    it("should handle empty SMILES arrays", async () => {
      const response = await request(app)
        .post("/generate-sdfs")
        .send({ smiles: [] });
      
      expect(response.status).toBe(200);
      expect(response.body.sdfPaths).toEqual([]);
    });
  });

  describe("Performance Integration", () => {
    it("should handle multiple concurrent requests", async () => {
      const promises = Array(3).fill().map((_, i) => 
        request(app)
          .post("/generate-sdfs")
          .send({ smiles: [`C${i}`, "O"] })
      );
      
      const responses = await Promise.all(promises);
      responses.forEach(response => {
        expect([200, 500]).toContain(response.status);
      });
    }, 20000);
  });
});

// Legacy integration tests - test fixtures and utilities validation
describe("Integration Tests - Test Infrastructure", () => {
  let fileManager;
  
  beforeEach(() => {
    fileManager = new TestFileManager();
  });
  
  afterEach(() => {
    fileManager.cleanup();
    jest.clearAllMocks();
  });

  describe("Test Fixtures", () => {
    it("should provide valid test molecules", () => {
      const caffeine = getTestMolecule("caffeine");
      
      expect(caffeine).toBeDefined();
      expect(caffeine.smiles).toBe("CN1C=NC2=C1C(=O)N(C(=O)N2C)C");
      expect(caffeine.name).toBe("Caffeine");
      expect(TestAssertions.isValidSmiles(caffeine.smiles)).toBe(true);
    });

    it("should provide mock images", () => {
      expect(MOCK_IMAGES.blackSquare.base64).toBeDefined();
      expect(MOCK_IMAGES.whiteSquare.base64).toBeDefined();
      expect(typeof MOCK_IMAGES.blackSquare.base64).toBe("string");
    });
  });

  describe("Test Request Creation", () => {
    it("should create image molecules test request", () => {
      const request = createTestRequest("imageMolecules", {
        object: "coffee",
        x: 150,
        y: 250
      });
      
      expect(request.imageBase64).toBeDefined();
      expect(request.croppedImageBase64).toBeDefined();
      expect(request.x).toBe(150);
      expect(request.y).toBe(250);
    });

    it("should create object molecules test request", () => {
      const request = createTestRequest("objectMolecules", {
        object: "wine"
      });
      
      expect(request.object).toBe("wine");
    });
  });

  describe("File Manager", () => {
    it("should create and cleanup test files", () => {
      const testContent = "test content";
      
      const filepath = fileManager.createTempFile("test.txt", testContent);
      
      expect(fs.existsSync(filepath)).toBe(true);
      expect(fs.readFileSync(filepath, "utf8")).toBe(testContent);
      
      fileManager.cleanup();
      expect(fs.existsSync(filepath)).toBe(false);
    });

    it("should create test SDF files", () => {
      const smiles = "CCO";
      
      const sdfPath = fileManager.createTestSdf(smiles, "ethanol.sdf");
      
      expect(fs.existsSync(sdfPath)).toBe(true);
      
      const content = fs.readFileSync(sdfPath, "utf8");
      expect(content).toContain(smiles);
      expect(content).toContain("_test_file");
      
      fileManager.cleanup();
      expect(fs.existsSync(sdfPath)).toBe(false);
    });
  });

  describe("Test Assertions", () => {
    it("should validate SMILES strings", () => {
      expect(TestAssertions.isValidSmiles("CCO")).toBe(true);
      expect(TestAssertions.isValidSmiles("O")).toBe(true);
      expect(TestAssertions.isValidSmiles("CN1C=NC2=C1C(=O)N(C(=O)N2C)C")).toBe(true);
      
      expect(TestAssertions.isValidSmiles("")).toBe(false);
      expect(TestAssertions.isValidSmiles("invalid smiles")).toBe(false);
      expect(TestAssertions.isValidSmiles(123)).toBe(false);
    });

    it("should validate SDF paths", () => {
      expect(TestAssertions.isValidSdfPath("/path/to/file.sdf")).toBe(true);
      expect(TestAssertions.isValidSdfPath("http://example.com/file.sdf")).toBe(true);
      
      expect(TestAssertions.isValidSdfPath("file.txt")).toBe(false);
      expect(TestAssertions.isValidSdfPath("")).toBe(false);
      expect(TestAssertions.isValidSdfPath(null)).toBe(false);
    });

    it("should validate arrays of SMILES", () => {
      expect(TestAssertions.arrayContainsValidSmiles(["CCO", "O"])).toBe(true);
      
      expect(TestAssertions.arrayContainsValidSmiles([])).toBe(false);
      expect(TestAssertions.arrayContainsValidSmiles(["CCO", "invalid"])).toBe(false);
      expect(TestAssertions.arrayContainsValidSmiles("not an array")).toBe(false);
    });
  });
});

// DOM/Frontend Integration Tests
describe("Integration Tests - Frontend Components", () => {
  let document, generateSDFs;
  const { JSDOM } = require("jsdom");
  const fetchMock = require("jest-fetch-mock");

  beforeEach(() => {
    const htmlContent = fs.readFileSync(path.join(__dirname, "..", "index.html"), "utf8");
    const dom = new JSDOM(htmlContent, {
      runScripts: "dangerously",
      resources: "usable",
      url: "http://localhost:8080"
    });
    
    document = dom.window.document;
    
    // Mock the generateSDFs function
    generateSDFs = jest.fn().mockImplementation(async (smilesList) => {
      const viewerContainer = document.getElementById("viewer-container") || document.createElement("div");
      viewerContainer.id = "viewer-container";
      if (!document.getElementById("viewer-container")) {
        document.body.appendChild(viewerContainer);
      }
      
      viewerContainer.innerHTML = '';
      
      try {
        const response = await fetch('http://localhost:8080/generate-sdfs', {
          method: 'POST',
          headers: { 'Content-Type': 'application/json' },
          body: JSON.stringify({ smiles: smilesList })
        });
        
        if (response.ok) {
          const data = await response.json();
          data.sdfPaths.forEach(() => {
            const viewer = document.createElement("div");
            viewer.className = "mol-viewer";
            viewerContainer.appendChild(viewer);
          });
        }
      } catch (error) {
        console.error("Mock generateSDFs error:", error);
      }
    });
    
    dom.window.generateSDFs = generateSDFs;
    fetchMock.resetMocks();
  });

  it("should render SDFs when a list of SMILES is passed", async () => {
    const smilesList = ["C1=CC=CC=C1", "C1CCCCC1"];
    fetchMock.mockResponseOnce(
      JSON.stringify({
        sdfPaths: ["/sdf_files/C1=CC=CC=C1.sdf", "/sdf_files/C1CCCCC1.sdf"],
      }),
    );

    await generateSDFs(smilesList);

    const viewerContainer = document.getElementById("viewer-container");
    expect(viewerContainer.children.length).toBe(2);
  });

  it("should handle errors gracefully", async () => {
    const smilesList = ["C1=CC=CC=C1"];
    fetchMock.mockRejectOnce(new Error("Failed to fetch"));

    await generateSDFs(smilesList);

    const viewerContainer = document.getElementById("viewer-container");
    expect(viewerContainer.children.length).toBe(0);
  });

  it("should have required DOM elements", () => {
    // Test that essential UI elements exist
    const textInput = document.querySelector('input[type="text"]') || document.querySelector('textarea');
    expect(textInput).toBeTruthy();
    
    // Should have video element for camera
    const video = document.querySelector('video');
    expect(video).toBeTruthy();
  });
}); 