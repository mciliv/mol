// server.js - Clean modular architecture
  process.env.NODE_ENV ||= 'development';

// ==================== IMPORTS ====================
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const HttpsServer = require("./https-server");
const AIAnalyzer = require("./ai-analyzer");
const MolecularProcessor = require("./molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("./schemas");

// ==================== CONFIGURATION ====================
const app = express();
const PORT = process.env.PORT || 8080;
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;

// Initialize modules
const aiAnalyzer = new AIAnalyzer(process.env.OPENAI_API_KEY);
const molecularProcessor = new MolecularProcessor();

// ==================== MIDDLEWARE ====================
app.use(cors());
app.use(express.json({ limit: '50mb' }));
app.use(express.static(__dirname));
app.use('/sdf_files', express.static(path.join(__dirname, 'sdf_files')));

// ==================== ROUTES ====================

// Image analysis route
app.post("/image-molecules", async (req, res) => {
  try {
    const { imageBase64, croppedImageBase64, x, y } = req.body;
    
    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }

    const result = await aiAnalyzer.analyzeImage(imageBase64, croppedImageBase64, x, y);
    res.json({ output: result });
  } catch (error) {
    console.error("Image analysis error:", error);
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route
app.post("/object-molecules", async (req, res) => {
  try {
    const { object } = req.body;
    
    if (!object) {
      return res.status(400).json({ error: "No object description provided" });
    }

    const result = await aiAnalyzer.analyzeText(object);
    res.json({ output: result });
  } catch (error) {
    console.error("Text analysis error:", error);
    res.status(500).json({ error: error.message });
  }
});

// SDF generation route
app.post('/generate-sdfs', async (req, res) => {
  try {
    const { smiles, overwrite = false } = req.body;
    
    if (!smiles || !Array.isArray(smiles)) {
      return res.status(400).json({ error: "smiles array is required" });
    }

    const result = await molecularProcessor.processSmiles(smiles, overwrite);
    
    res.json({
      sdfPaths: result.sdfPaths,
      errors: result.errors,
      skipped: result.skipped,
      message: `Generated ${result.sdfPaths.length} 3D structures from ${smiles.length} SMILES`
    });
  } catch (error) {
    console.error("SDF generation error:", error);
    res.status(500).json({ error: error.message });
  }
});

// ==================== DEVELOPMENT MIDDLEWARE ====================
if (process.env.NODE_ENV === "development") {
  const livereload = require("livereload");
  const connectLivereload = require("connect-livereload");

  try {
    const liveReloadServer = livereload.createServer({
      exts: ['html', 'css', 'js'],
      ignore: ['node_modules/**', 'sdf_files/**', '*.log'],
      port: 35732  // Use different port to avoid conflicts
    });
    liveReloadServer.watch(__dirname);

    app.use(connectLivereload());

    liveReloadServer.server.once("connection", () => {
      setTimeout(() => {
        liveReloadServer.refresh("/");
      }, 100);
    });
    
    console.log("✅ LiveReload server started");
  } catch (err) {
    console.log("⚠️ LiveReload server failed to start:", err.message);
    console.log("🔄 Continuing without LiveReload...");
  }
}

// Static routes
app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "index.html"));
});

// Request logging middleware
app.use((req, res, next) => {
  console.log(`Incoming request: ${req.method} ${req.url}`);
  next();
});

// ==================== SERVER STARTUP ====================
const isCloudFunction = process.env.FUNCTION_NAME || process.env.FUNCTION_TARGET || process.env.K_SERVICE || process.env.GOOGLE_CLOUD_PROJECT || process.env.GCP_PROJECT;
const isNetlify = process.env.NETLIFY;
const isTestMode = process.env.NODE_ENV === 'test' || process.env.JEST_WORKER_ID;
const isServerless = isCloudFunction || isNetlify;

if (!isServerless && !isTestMode) {
  // Local development mode
  app.listen(PORT, '0.0.0.0', () => {
    console.log(`🚀 HTTP server running on http://0.0.0.0:${PORT}`);
  });

  // Start HTTPS server for development
  if (process.env.NODE_ENV !== 'production') {
    const httpsServer = new HttpsServer(app, HTTPS_PORT);
    httpsServer.start();
  }
} else {
  // Serverless mode
  if (isCloudFunction) {
    console.log(`☁️ Running in Cloud Functions mode`);
  } else if (isNetlify) {
    console.log(`◈ Running in Netlify mode`);
  }
  
  if (process.env.NODE_ENV === 'production' && !process.env.OPENAI_API_KEY) {
    console.error('❌ OPENAI_API_KEY environment variable is required for production');
    process.exit(1);
  }
}

// Always export the app for Cloud Functions
module.exports = app;