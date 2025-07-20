// server.js - Clean modular architecture
process.env.NODE_ENV ||= "development";

// ==================== IMPORTS ====================
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const HttpsServer = require("./https-server");
const AtomPredictor = require("./AtomPredictor");
const MolecularProcessor = require("./molecular-processor");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
} = require("./schemas");

// ==================== CONFIGURATION ====================
const app = express();
const DEFAULT_PORT = 8080;
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;

// Utility to attempt cleanup of processes using our ports
const attemptPortCleanup = async (port) => {
  try {
    const { execSync } = require("child_process");
    
    // Try to find and kill processes using the port (Unix/macOS)
    if (process.platform !== "win32") {
      try {
        const result = execSync(`lsof -ti:${port}`, { encoding: "utf8", stdio: "pipe" });
        const pids = result.trim().split("\n").filter(Boolean);
        
        if (pids.length > 0) {
          console.log(`🧹 Found ${pids.length} process(es) using port ${port}, attempting cleanup...`);
          
          for (const pid of pids) {
            try {
              execSync(`kill -9 ${pid}`, { stdio: "pipe" });
              console.log(`   ✅ Killed process ${pid}`);
            } catch (e) {
              console.log(`   ⚠️ Could not kill process ${pid} (may not have permission)`);
            }
          }
          
          // Wait a moment for cleanup
          await new Promise(resolve => setTimeout(resolve, 500));
          return true;
        }
      } catch (e) {
        // No processes found or lsof not available
        return false;
      }
    }
    return false;
  } catch (error) {
    console.log(`⚠️ Port cleanup failed: ${error.message}`);
    return false;
  }
};

const findAvailablePort = async (startPort) => {
  const net = require("net");

  const isPortAvailable = (port) => {
    return new Promise((resolve) => {
      const tester = net
        .createServer()
        .once("error", () => resolve(false))
        .once("listening", () => {
          tester.once("close", () => resolve(true)).close();
        })
        .listen(port);
    });
  };

  let port = startPort;
  while (port < startPort + 100) {
    // Try up to 100 ports
    if (await isPortAvailable(port)) {
      return port;
    }
    port++;
  }
  throw new Error(
    `No available ports found between ${startPort} and ${startPort + 100}`,
  );
};

const PORT = process.env.PORT || DEFAULT_PORT;

// Initialize modules
const atomPredictor = new AtomPredictor(process.env.OPENAI_API_KEY);
const molecularProcessor = new MolecularProcessor();

// ==================== MIDDLEWARE ====================
app.use(cors());
app.use(express.json({ limit: "50mb" }));

// In-memory user storage (for demo - replace with database in production)
const users = new Map(); // deviceToken -> user data

// ==================== DEVELOPMENT MIDDLEWARE ====================
// Live reload disabled - using external tools if needed
if (false && process.env.NODE_ENV === "development") {
  const livereload = require("livereload");
  const connectLivereload = require("connect-livereload");
  const net = require("net");

  const isPortAvailable = (port) => {
    return new Promise((resolve) => {
      const tester = net
        .createServer()
        .once("error", () => resolve(false))
        .once("listening", () => {
          tester.once("close", () => resolve(true)).close();
        })
        .listen(port);
    });
  };

  const startLiveReload = async (port) => {
    try {
      const available = await isPortAvailable(port);
      if (!available) {
        console.log(`LiveReload port ${port} is already in use`);
        if (port === 35729) {
          console.log("Trying alternative port 35730...");
          return startLiveReload(35730);
        } else {
          console.log("Continuing without LiveReload...");
          return;
        }
      }

      const liveReloadServer = livereload.createServer({
        exts: ["html", "css", "js"],
        ignore: ["node_modules/**", "sdf_files/**", "*.log"],
        port: port,
      });

      // Add error handler for unexpected errors
      liveReloadServer.server.on("error", (err) => {
        console.error("LiveReload server error:", err.message);
      });

      // Only proceed if server starts successfully
      liveReloadServer.server.once("listening", () => {
        liveReloadServer.watch(__dirname);
        app.use(connectLivereload());

        liveReloadServer.server.once("connection", () => {
          setTimeout(() => {
            liveReloadServer.refresh("/");
          }, 100);
        });

        console.log(`LiveReload server started on port ${port}`);
      });
    } catch (err) {
      console.log("LiveReload server failed to start:", err.message);
      console.log("Continuing without LiveReload...");
    }
  };

  // Start LiveReload with fallback port
  startLiveReload(35729);
}

app.use(express.static(__dirname));
app.use("/sdf_files", express.static(path.join(__dirname, "sdf_files")));

// ==================== PAYMENT ROUTES ====================

// Stripe configuration endpoint
app.get("/stripe-config", (req, res) => {
  res.json({
    publishableKey: process.env.STRIPE_PUBLISHABLE_KEY || 'pk_test_demo_key_for_development'
  });
});

// Setup payment method endpoint
app.post("/setup-payment-method", async (req, res) => {
  try {
    const { payment_method, device_info, name } = req.body;
    
    if (!payment_method || !device_info) {
      return res.status(400).json({ error: "Payment method and device info required" });
    }

    // Generate device token
    const deviceToken = Buffer.from(`${device_info}-${Date.now()}-${Math.random()}`).toString('base64').replace(/[+/=]/g, '');
    
    // Store user info (in production, this would go to a database)
    const userData = {
      deviceToken,
      paymentMethodId: payment_method,
      deviceInfo: device_info,
      name: name || null,
      usage: 0,
      createdAt: new Date(),
      lastUsed: new Date()
    };
    
    users.set(deviceToken, userData);
    
    // In production, you would:
    // 1. Create customer in Stripe
    // 2. Attach payment method to customer
    // 3. Handle 3D Secure if needed
    // For demo, we'll just return success
    
    console.log(`✅ New user setup: ${name || 'Anonymous'} with device ${deviceToken.substring(0, 8)}...`);
    
    res.json({
      success: true,
      device_token: deviceToken,
      requires_action: false // Set to true if 3D Secure needed
    });
    
  } catch (error) {
    console.error("Payment setup error:", error);
    res.status(500).json({ error: "Failed to setup payment method" });
  }
});

// Validate payment method endpoint
app.post("/validate-payment", (req, res) => {
  try {
    const { device_token } = req.body;
    
    if (!device_token) {
      return res.status(400).json({ error: "Device token required" });
    }
    
    const user = users.get(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    // Update last used
    user.lastUsed = new Date();
    
    res.json({
      valid: true,
      user: {
        name: user.name,
        usage: user.usage,
        device_token: user.deviceToken
      }
    });
    
  } catch (error) {
    console.error("Payment validation error:", error);
    res.status(500).json({ error: "Failed to validate payment" });
  }
});

// Increment usage endpoint
app.post("/increment-usage", (req, res) => {
  try {
    const { device_token } = req.body;
    
    if (!device_token) {
      return res.status(400).json({ error: "Device token required" });
    }
    
    const user = users.get(device_token);
    if (!user) {
      return res.status(404).json({ error: "User not found" });
    }
    
    user.usage++;
    user.lastUsed = new Date();
    
    console.log(`📊 Usage: ${user.name || 'Anonymous'} - ${user.usage} analyses`);
    
    res.json({
      usage: user.usage
    });
    
  } catch (error) {
    console.error("Usage increment error:", error);
    res.status(500).json({ error: "Failed to increment usage" });
  }
});

// ==================== ANALYSIS ROUTES ====================

// Image analysis route
app.post("/image-molecules", async (req, res) => {
  try {
    // Validate input schema
    const validation = ImageMoleculeSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
      });
    }

    const {
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    } = req.body;

    if (!imageBase64) {
      return res.status(400).json({ error: "No image data provided" });
    }

    const result = await atomPredictor.analyzeImage(
      imageBase64,
      croppedImageBase64,
      x,
      y,
      cropMiddleX,
      cropMiddleY,
      cropSize,
    );
    res.json({ output: result });
  } catch (error) {
    console.error("Image analysis error:", error);

    // Provide more specific error messages
    let errorMessage = error.message;
    let statusCode = 500;

    if (error.message.includes("network") || error.message.includes("fetch")) {
      errorMessage =
        "Network error: Unable to connect to AI service. Please check your internet connection.";
      statusCode = 503;
    } else if (
      error.message.includes("API key") ||
      error.message.includes("authentication")
    ) {
      errorMessage = "Authentication error: Invalid or missing API key.";
      statusCode = 401;
    } else if (
      error.message.includes("rate limit") ||
      error.message.includes("quota")
    ) {
      errorMessage = "Rate limit exceeded: Please try again later.";
      statusCode = 429;
    } else if (error.message.includes("timeout")) {
      errorMessage =
        "Request timeout: The AI service is taking too long to respond.";
      statusCode = 408;
    }

    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
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
        details: validation.error.issues,
      });
    }

    const { object } = req.body;

    if (!object) {
      return res.status(400).json({ error: "No object description provided" });
    }

    const result = await atomPredictor.analyzeText(object);
    res.json({ output: result });
  } catch (error) {
    console.error("Text analysis error:", error);

    // Provide more specific error messages
    let errorMessage = error.message;
    let statusCode = 500;

    if (error.message.includes("network") || error.message.includes("fetch")) {
      errorMessage =
        "Network error: Unable to connect to AI service. Please check your internet connection.";
      statusCode = 503;
    } else if (
      error.message.includes("API key") ||
      error.message.includes("authentication")
    ) {
      errorMessage = "Authentication error: Invalid or missing API key.";
      statusCode = 401;
    } else if (
      error.message.includes("rate limit") ||
      error.message.includes("quota")
    ) {
      errorMessage = "Rate limit exceeded: Please try again later.";
      statusCode = 429;
    } else if (error.message.includes("timeout")) {
      errorMessage =
        "Request timeout: The AI service is taking too long to respond.";
      statusCode = 408;
    }

    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
  }
});

// SDF generation route
app.post("/generate-sdfs", async (req, res) => {
  try {
    // Validate input schema
    const validation = SdfGenerationSchema.safeParse(req.body);
    if (!validation.success) {
      return res.status(400).json({
        error: "Invalid input data",
        details: validation.error.issues,
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
      message: `Generated ${result.sdfPaths.length} 3D structures from ${smiles.length} SMILES`,
    });
  } catch (error) {
    console.error("SDF generation error:", error);

    // Provide more specific error messages
    let errorMessage = error.message;
    let statusCode = 500;

    if (error.message.includes("network") || error.message.includes("fetch")) {
      errorMessage =
        "Network error: Unable to connect to molecular service. Please check your internet connection.";
      statusCode = 503;
    } else if (
      error.message.includes("file system") ||
      error.message.includes("permission")
    ) {
      errorMessage =
        "File system error: Unable to create SDF files. Please check directory permissions.";
      statusCode = 500;
    } else if (error.message.includes("timeout")) {
      errorMessage =
        "Request timeout: Molecular structure generation is taking too long.";
      statusCode = 408;
    }

    res.status(statusCode).json({
      error: errorMessage,
      details: process.env.NODE_ENV === "development" ? error.stack : undefined,
    });
  }
});

// Static routes
app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "index.html"));
});

// Request logging middleware
app.use((req, res, next) => {
  // Skip logging Chrome DevTools discovery requests
  if (!req.url.includes(".well-known/appspecific/com.chrome.devtools")) {
    console.log(`Incoming request: ${req.method} ${req.url}`);
  }

  // Handle Chrome DevTools discovery request
  if (req.url === "/.well-known/appspecific/com.chrome.devtools.json") {
    return res.status(404).json({ error: "Not found" });
  }

  next();
});

// ==================== SERVER STARTUP ====================
const isCloudFunction =
  process.env.FUNCTION_NAME ||
  process.env.FUNCTION_TARGET ||
  process.env.K_SERVICE ||
  process.env.GOOGLE_CLOUD_PROJECT ||
  process.env.GCP_PROJECT;
const isNetlify = process.env.NETLIFY;
const isTestMode =
  process.env.NODE_ENV === "test" || process.env.JEST_WORKER_ID;
const isServerless = isCloudFunction || isNetlify;

// Store server instances for cleanup
let httpServer;
let httpsServerInstance;

if (!isServerless && !isTestMode) {
  // Local development mode
  const startServer = async () => {
    try {
      let actualPort = PORT;

      // If default port is in use, try cleanup and then find an available port
      if (PORT === DEFAULT_PORT) {
        try {
          // First check if port is available
          if (!(await (async () => {
            const net = require("net");
            return new Promise((resolve) => {
              const server = net.createServer();
              server.listen(PORT, "0.0.0.0", () => {
                server.once("close", () => resolve(true));
                server.close();
              });
              server.on("error", () => resolve(false));
            });
          })())) {
            console.log(`⚠️ Port ${PORT} is in use, attempting cleanup...`);
            const cleanupSuccessful = await attemptPortCleanup(PORT);
            
            if (cleanupSuccessful) {
              console.log(`✅ Port ${PORT} cleanup completed, retrying...`);
              // Give a moment for the port to be fully released
              await new Promise(resolve => setTimeout(resolve, 1000));
            }
          }
          
          actualPort = await findAvailablePort(PORT);
          if (actualPort !== PORT) {
            console.log(
              `⚠️  Port ${PORT} is in use, using port ${actualPort} instead`,
            );
          }
        } catch (error) {
          console.error(`❌ Could not find available port: ${error.message}`);
          console.log(`💡 Try: pkill -f "node.*server.js" or use a different port`);
          process.exit(1);
        }
      }

      httpServer = app.listen(actualPort, "0.0.0.0", () => {
        console.log(`✅ HTTP server running on http://localhost:${actualPort}`);
        console.log(
          `📱 Mobile access: http://$(ifconfig | grep 'inet ' | grep -v 127.0.0.1 | awk '{print $2}' | head -1):${actualPort}`,
        );
      });
    } catch (error) {
      console.error(`❌ Failed to start server: ${error.message}`);
      process.exit(1);
    }
  };

  startServer()
    .then(() => {
      // Handle port binding errors
      httpServer.on("error", (error) => {
        if (error.code === "EADDRINUSE") {
          console.error(`❌ Port ${PORT} is already in use`);
          console.log(`💡 Solutions:`);
          console.log(
            `   1. Kill existing process: pkill -f "node.*server.js"`,
          );
          console.log(`   2. Use different port: PORT=8081 npm start`);
          console.log(`   3. Check what's using the port: lsof -i :${PORT}`);
          process.exit(1);
        } else if (error.code === "EACCES") {
          console.error(`❌ Permission denied: Cannot bind to port ${PORT}`);
          console.log(`💡 Try using a port > 1024 or run with sudo`);
          process.exit(1);
        } else {
          console.error("❌ Server error:", error.message);
          console.log(`💡 Check your network configuration and try again`);
          process.exit(1);
        }
      });
    })
    .catch((error) => {
      console.error(`❌ Failed to start server: ${error.message}`);
      process.exit(1);
    });

  // Start HTTPS server for development
  if (process.env.NODE_ENV !== "production") {
    const startHttpsServer = async () => {
      try {
        const httpsServer = new HttpsServer(app, HTTPS_PORT);
        httpsServerInstance = await httpsServer.start();
        
        if (httpsServerInstance) {
          console.log("✅ HTTPS server started successfully");
          
          // Handle HTTPS server errors after startup
          httpsServerInstance.on("error", (error) => {
            console.error("❌ HTTPS server error after startup:", error.message);
            console.log("💡 HTTPS server will continue running if possible");
          });
        } else {
          console.log("⚠️ HTTPS server not started - continuing with HTTP only");
        }
      } catch (error) {
        console.error("❌ Failed to start HTTPS server:", error.message);
        console.log("💡 Continuing with HTTP server only");
      }
    };

    // Start HTTPS server after a short delay to avoid port conflicts
    setTimeout(startHttpsServer, 1000);
  }
} else {
  // Serverless mode
  if (isCloudFunction) {
    console.log(`Running in Cloud Functions mode`);
  } else if (isNetlify) {
    console.log(`Running in Netlify mode`);
  }

  if (process.env.NODE_ENV === "production" && !process.env.OPENAI_API_KEY) {
    console.error(
      "OPENAI_API_KEY environment variable is required for production",
    );
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
            console.log(`${name} server closed`);
            resolve();
          });
        } else {
          resolve();
        }
      });
    };

    Promise.all([
      closeServer(httpServer, "HTTP"),
      closeServer(httpsServerInstance, "HTTPS"),
    ]).then(() => {
      console.log("Shutdown complete");
      process.exit(0);
    });

    // Force exit after 5 seconds
    setTimeout(() => {
      console.error("Forced shutdown after timeout");
      process.exit(1);
    }, 5000);
  };

  // Listen for termination signals
  process.on("SIGTERM", () => gracefulShutdown("SIGTERM"));
  process.on("SIGINT", () => gracefulShutdown("SIGINT"));
  process.on("SIGUSR2", () => gracefulShutdown("SIGUSR2")); // Nodemon uses this
}

// Global error handlers
process.on("unhandledRejection", (reason, promise) => {
  console.error("❌ Unhandled Rejection at:", promise, "reason:", reason);
  console.log(
    "💡 This usually indicates a network or API error. Check your internet connection and API keys.",
  );
});

process.on("uncaughtException", (error) => {
  console.error("❌ Uncaught Exception:", error.message);
  console.log("💡 Application crashed. Check the error details above.");
  process.exit(1);
});

// Always export the app for Cloud Functions
module.exports = app;
