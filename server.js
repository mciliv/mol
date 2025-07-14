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
    // Validate input schema
    const validation = ImageMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({ 
        error: "Invalid input data", 
        details: validation.error.issues 
      });
    }

    const { imageBase64, croppedImageBase64, x, y, cropMiddleX, cropMiddleY, cropSize } = req.body;
    
    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }

    const result = await aiAnalyzer.analyzeImage(imageBase64, croppedImageBase64, x, y, cropMiddleX, cropMiddleY, cropSize);
    res.json({ output: result });
  } catch (error) {
    console.error("Image analysis error:", error);
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route
app.post("/object-molecules", async (req, res) => {
  try {
    // Validate input schema
    const validation = TextMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({ 
        error: "Invalid input data", 
        details: validation.error.issues 
      });
    }

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
    // Validate input schema
    const validation = SdfGenerationSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({ 
        error: "Invalid input data", 
        details: validation.error.issues 
      });
    }

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
      port: 0  // Let system assign an available port automatically
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

// Store server instances for cleanup
let httpServer;
let httpsServerInstance;

if (!isServerless && !isTestMode) {
  // Local development mode
  httpServer = app.listen(PORT, '0.0.0.0', () => {
    console.log(`🚀 HTTP server running on http://0.0.0.0:${PORT}`);
  });
  
  // Handle port binding errors
  httpServer.on('error', (error) => {
    if (error.code === 'EADDRINUSE') {
      console.error(`❌ Port ${PORT} is already in use`);
      console.log(`💡 Try: pkill -f "node.*server.js" or use a different port`);
      process.exit(1);
    } else {
      console.error('Server error:', error);
      process.exit(1);
    }
  });

  // Start HTTPS server for development
  if (process.env.NODE_ENV !== 'production') {
    const httpsServer = new HttpsServer(app, HTTPS_PORT);
    httpsServerInstance = httpsServer.start();
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

// Graceful shutdown handling for nodemon restarts
if (!isServerless && !isTestMode) {
  const gracefulShutdown = (signal) => {
    console.log(`\n${signal} received: closing HTTP/HTTPS servers gracefully`);
    
    const closeServer = (server, name) => {
      return new Promise((resolve) => {
        if (server) {
          server.close(() => {
            console.log(`✅ ${name} server closed`);
            resolve();
          });
        } else {
          resolve();
        }
      });
    };
    
    Promise.all([
      closeServer(httpServer, 'HTTP'),
      closeServer(httpsServerInstance, 'HTTPS')
    ]).then(() => {
      console.log('👋 Shutdown complete');
      process.exit(0);
    });
    
    // Force exit after 5 seconds
    setTimeout(() => {
      console.error('⚠️ Forced shutdown after timeout');
      process.exit(1);
    }, 5000);
  };
  
  // Listen for termination signals
  process.on('SIGTERM', () => gracefulShutdown('SIGTERM'));
  process.on('SIGINT', () => gracefulShutdown('SIGINT'));
  process.on('SIGUSR2', () => gracefulShutdown('SIGUSR2')); // Nodemon uses this
}

// Always export the app for Cloud Functions
module.exports = app;