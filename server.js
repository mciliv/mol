// server.js (with dynamic LAN IP detection for HTTPS)
process.env.NODE_ENV ||= 'development';

// Load environment variables
if (process.env.NODE_ENV !== 'production') {
  try {
    require('dotenv').config();
  } catch (err) {
    console.log('⚠️ dotenv not found, using environment variables directly');
  }
}

// ==================== IMPORTS ====================
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const { spawn, execSync } = require("child_process");
const https = require("https");
const os = require("os");
const OpenAI = require("openai");
const { zodTextFormat } = require("openai/helpers/zod");
const {
  ObjectIdentificationSchema,
  ImageMoleculeSchema,
  TextMoleculeSchema,
  ListMoleculesTextRequestSchema,
} = require("./schemas");

// ==================== CONFIGURATION ====================
const app = express();
const PORT = process.env.PORT || 8080;  // Cloud Functions default, was 3000
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;  // Configurable HTTPS port
const SDF_DIR = path.join(__dirname, "sdf_files");

// Ensure SDF directory exists
if (!fs.existsSync(SDF_DIR)) fs.mkdirSync(SDF_DIR, { recursive: true });

// ==================== UTILITY FUNCTIONS ====================
function getLocalIPAddress() {
  const interfaces = os.networkInterfaces();
  for (const name of Object.keys(interfaces)) {
    for (const iface of interfaces[name]) {
      if (iface.family === 'IPv4' && !iface.internal) {
        return iface.address;
      }
    }
  }
  return '127.0.0.1';
}

function generateSelfSignedCert(ip) {
  const certDir = path.join(__dirname, 'certs');
  const keyPath = path.join(certDir, 'key.pem');
  const certPath = path.join(certDir, 'cert.pem');

  if (fs.existsSync(keyPath) && fs.existsSync(certPath)) {
    console.log("✅ Using existing SSL certificates");
    return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
  }

  if (!fs.existsSync(certDir)) fs.mkdirSync(certDir, { recursive: true });

  console.log("🔐 Generating self-signed SSL certificates for development...");

  try {
    execSync(`openssl genrsa -out ${keyPath} 2048`, { stdio: 'inherit' });
    const opensslCmd = [
      'req', '-new', '-x509', '-key', keyPath, '-out', certPath, '-days', '365',
      '-subj', '/C=US/ST=Dev/L=Local/O=MolecularReality/CN=localhost',
      '-addext', `subjectAltName=DNS:localhost,IP:127.0.0.1,IP:0.0.0.0,IP:${ip}`
    ];
    execSync(`openssl ${opensslCmd.join(' ')}`, { stdio: 'inherit' });

    console.log("✅ SSL certificates generated successfully!");
    return { key: fs.readFileSync(keyPath), cert: fs.readFileSync(certPath) };
  } catch (err) {
    console.error("❌ Failed to generate SSL certificates:", err.message);
    return null;
  }
}

// ==================== AI/CHEMICAL ANALYSIS FUNCTIONS ====================
function sdf(s, overwrite) {
  let command = "python";
  let args = ["sdf.py", s, "--dir", SDF_DIR];
  if (overwrite) args.push("--overwrite");
  return { command, args };
}

function crystal(s, overwrite) {
  let command = "python";
  let args = ["crystal.py", s, "--dir", SDF_DIR];
  if (overwrite) args.push("--overwrite");
  return { command, args };
}

// Known mineral formulas that can be converted to crystal structures
const KNOWN_MINERALS = [
  'CaCO3', 'SiO2', 'Al2O3', 'FeS2', 'NaCl',
  'CaCO₃', 'SiO₂', 'Al₂O₃', 'FeS₂',
  'quartz', 'calcite', 'corundum', 'pyrite', 'halite', 'salt'
];

function isMineralFormula(chemical) {
  // Check if it's a known mineral
  if (KNOWN_MINERALS.some(mineral => 
    chemical.toLowerCase() === mineral.toLowerCase())) {
    return true;
  }
  
  // Exclude pure organic molecular formulas (CxHy, CxHyN, etc.) - these should be SMILES
  if (/^C\d*H\d*(N\d*|O\d*|S\d*)*$/.test(chemical)) {
    return false;
  }
  
  // Check if it looks like a mineral formula (contains metals/metalloids)
  const hasMetals = /[A-Z][a-z]?/.test(chemical) && 
    !/^[CHONPS]+\d*$/.test(chemical); // Not just organic elements
  const hasNoSmilesChars = !/[\[\]()=#+\-\.@:\/\\%]/.test(chemical);
  
  return hasMetals && hasNoSmilesChars;
}

function getChemicalProcessor(chemical) {
  // Always try SMILES first, crystal as backup
  return { type: 'smiles', processor: sdf, backup: crystal };
}

function findGeneratedSdfFile(chemical, sdfDir) {
  // Try different possible filenames that the Python scripts might create
  const possibleFilenames = [
    `${chemical}.sdf`,  // Original filename
    `${chemical.replace(/[^a-zA-Z0-9]/g, '_')}.sdf`,  // Safe filename
  ];
  
  // For minerals, also try normalized names
  if (isMineralFormula(chemical)) {
    const normalized = parseMineral(chemical);
    if (normalized !== chemical) {
      possibleFilenames.push(`${normalized}.sdf`);
    }
  }
  
  // Check which file actually exists
  for (const filename of possibleFilenames) {
    const fullPath = path.join(sdfDir, filename);
    if (fs.existsSync(fullPath)) {
      return `/sdf_files/${filename}`;  // Return web path
    }
  }
  
  return null;
}

function parseMineral(formula) {
  // Simple mineral name normalization (matches crystal.py logic)
  const variations = {
    'quartz': 'SiO2',
    'calcite': 'CaCO3', 
    'corundum': 'Al2O3',
    'pyrite': 'FeS2',
    'halite': 'NaCl',
    'salt': 'NaCl',
    'CaCO₃': 'CaCO3',
    'SiO₂': 'SiO2', 
    'Al₂O₃': 'Al2O3',
    'FeS₂': 'FeS2'
  };
  
  return variations[formula.toLowerCase()] || formula;
}

// Configurable analysis modes
const ANALYSIS_MODES = {
  ORGANIC_ONLY: {
    name: "organic_only",
    instructions: `IMPORTANT: Return only valid SMILES notation (e.g., CCO for ethanol, C1=CC=CC=C1 for benzene).
Do NOT return molecular formulas (like C6H10O5, CaCO3, SiO2).
Do NOT return complex mineral formulas (like Mg3Si4O10(OH)2).
If the material contains inorganic compounds that cannot be represented in SMILES, focus on organic components only.
Each string must be parseable by RDKit/OpenBabel.`
  },
  
  MIXED_COMPOUNDS: {
    name: "mixed_compounds", 
    instructions: `Return chemical structures in the most appropriate format:
- For organic molecules: use SMILES notation (e.g., CCO for ethanol, C1=CC=CC=C1 for benzene)
- For simple inorganics that have SMILES: use SMILES (e.g., O for water, [Na+].[Cl-] for salt)
- For minerals/crystals: use molecular formulas (e.g., CaCO3, SiO2, Mg3Si4O10(OH)2)
- For metals: use symbols (e.g., Au, Fe, Cu)
Each string should be chemically valid and parseable when possible.`
  },

  MINERALS_FOCUS: {
    name: "minerals_focus",
    instructions: `Focus on mineral and crystalline structures:
- Return molecular formulas for minerals (e.g., CaCO3, SiO2, Al2O3)
- Include crystal structure information when relevant
- For organic components in minerals: use simple molecular formulas
- Prioritize geological and mineralogical accuracy over SMILES compatibility.`
  }
};

// Current analysis mode - can be changed via environment variable or API
let CURRENT_MODE = process.env.ANALYSIS_MODE || 'MIXED_COMPOUNDS';

// Direct constant for current instructions - updates when mode changes
let SMILES_INSTRUCTIONS = (ANALYSIS_MODES[CURRENT_MODE] || ANALYSIS_MODES.MIXED_COMPOUNDS).instructions;

async function getSmilesForObject(object, imageBase64 = null, croppedImageBase64 = null) {
  const client = new OpenAI({ apiKey: process.env["OPENAI_API_KEY"] });
  const smilesRules = SMILES_INSTRUCTIONS;
  
  let text = `Identify: ${object}. List relevant chemical structures as VALID SMILES strings only.

${smilesRules}

Consider the specific material, compound, or substance shown.`;
  
  let content = [{ text: text, type: "input_text" }];
  
  // Add images if provided for better context
  if (imageBase64) {
    console.log(`🖼️  Using image context for enhanced SMILES analysis`);
    content.push({ detail: "high", type: "input_image", image_url: `data:image/jpeg;base64,${imageBase64}` });
    text = `Analyze the object "${object}" in this image. List relevant chemical structures as VALID SMILES strings only.

${smilesRules}

Use visual details to be more specific about materials, compounds, or substances shown.`;
    content[0].text = text;
  }
  
  if (croppedImageBase64) {
    console.log(`🔍 Using cropped region for focused chemical analysis`);
    content.push({ detail: "high", type: "input_image", image_url: `data:image/jpeg;base64,${croppedImageBase64}` });
  }

  const parsed = await client.responses.parse({
    input: [{ content: content, role: "user" }],
    model: "gpt-4o",
    text: { format: zodTextFormat(TextMoleculeSchema, "text_molecule_schema") },
  });
  
  // Handle cases where OpenAI returns null or no smiles
  if (!parsed.output_parsed || !parsed.output_parsed.smiles) {
    console.log("⚠️ No SMILES found in OpenAI response");
    return [];
  }
  
  return parsed.output_parsed.smiles;
}

async function identifyObjectFromImage(imageBase64, croppedImageBase64, x, y) {
  const client = new OpenAI({ apiKey: process.env["OPENAI_API_KEY"] });
  
  const identificationText = `You are a chemical analysis expert. The user clicked at coordinate (X: ${Math.round(x)}, Y: ${Math.round(y)}) in large image. Identify the specific substance, material, or chemical compound at that location. Be as specific as possible, and particularly in anyway that would effect its molecular structure - consider molecular composition, material type, active ingredients, or chemical structure when naming the object. Allowing human analysis st it'd contriute to human safety &  health`;

  // Save images for debugging
  fs.writeFileSync("image.jpg", imageBase64, "base64");
  fs.writeFileSync("cropped_image.jpg", croppedImageBase64, "base64");

  const objectIdentification = await client.responses.parse({
    input: [
      {
        content: [
          { text: identificationText, type: "input_text" },
          { detail: "high", type: "input_image", image_url: `data:image/jpeg;base64,${imageBase64}` },
          { detail: "high", type: "input_image", image_url: `data:image/jpeg;base64,${croppedImageBase64}` },
        ],
        role: "user",
      },
    ],
    model: "gpt-4o",
    text: { format: zodTextFormat(ObjectIdentificationSchema, "object_identification_schema") },
  });

  return objectIdentification.output_parsed.object;
}

// ==================== DEVELOPMENT MIDDLEWARE ====================
if (process.env.NODE_ENV === "development") {
  const livereload = require("livereload");
  const connectLivereload = require("connect-livereload");

  try {
    const liveReloadServer = livereload.createServer({
      exts: ['html', 'css', 'js'],
      ignore: ['node_modules/**', 'sdf_files/**', '*.log'],
      port: 35730  // Use different port to avoid conflicts
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

// ==================== MIDDLEWARE SETUP ====================
app.use((req, res, next) => {
  console.log(`Incoming request: ${req.method} ${req.url}`);
  next();
});

app.use(cors());
app.use(express.json({ limit: "50mb" }));
app.use(express.static(__dirname));
app.use("/sdf_files", express.static(SDF_DIR));
app.use("/favicon.ico", express.static(path.join(__dirname, "favicon.ico")));

// ==================== ROUTES ====================
// Static routes
app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "index.html"));
});

// Image analysis route
app.post("/image-molecules", async (req, res) => {
  const { imageBase64, croppedImageBase64, x, y } = req.body;

  if (!imageBase64) {
    return res.status(400).json({ error: "No image data provided" });
  }

  try {
    // Step 1: Identify the object from the image
    const identifiedObject = await identifyObjectFromImage(imageBase64, croppedImageBase64, x, y);
    console.log(`🔍 Identified object: ${identifiedObject}`);

    // Step 2: Get SMILES for the identified object (with image context)
    console.log(`🧪 Analyzing ${identifiedObject} for chemical structures...`);
    const smiles = await getSmilesForObject(identifiedObject, imageBase64, croppedImageBase64);

    const result = {
      object: identifiedObject,
      smiles: smiles
    };

    console.log(`✅ Final result: ${JSON.stringify(result, null, 2)}`);
    return res.json({ output: result });
  } catch (error) {
    console.error("Error:", error);
    if (error.name === "ZodError") {
      return res
        .status(400)
        .json({ error: "Invalid response format", details: error.errors });
    }
    res.status(500).json({ error: error.message });
  }
});

// Text analysis route
app.post("/object-molecules", async (req, res) => {
  try {
    const validatedData = ListMoleculesTextRequestSchema.parse(req.body);
    const { object } = validatedData;

    const smiles = await getSmilesForObject(object);
    return res.json({ output: { smiles } });
  } catch (error) {
    console.error("Error in /object-molecules handler:", error);
    if (error.name === "ZodError") {
      return res
        .status(400)
        .json({ error: "Invalid request data", details: error.errors });
    }
    return res.status(500).json({ error: error.message });
  }
});

// Analysis mode management routes
app.get("/analysis-mode", (req, res) => {
  const mode = ANALYSIS_MODES[CURRENT_MODE] || ANALYSIS_MODES.MIXED_COMPOUNDS;
  res.json({ 
    current_mode: CURRENT_MODE,
    available_modes: Object.keys(ANALYSIS_MODES),
    description: mode.name,
    instructions: mode.instructions
  });
});

app.post("/analysis-mode", (req, res) => {
  const { mode } = req.body;
  
  if (!mode || !ANALYSIS_MODES[mode]) {
    return res.status(400).json({ 
      error: "Invalid mode", 
      available_modes: Object.keys(ANALYSIS_MODES)
    });
  }
  
  // Update current mode (in production, this would be stored in database)
  CURRENT_MODE = mode;
  SMILES_INSTRUCTIONS = ANALYSIS_MODES[mode].instructions;
  process.env.ANALYSIS_MODE = mode;
  
  res.json({ 
    message: `Analysis mode changed to ${mode}`,
    current_mode: mode,
    instructions: ANALYSIS_MODES[mode].instructions
  });
});

// SDF generation route (updated to handle both SMILES and minerals)
app.post('/generate-sdfs', async (req, res) => {
  const { smiles, overwrite = false } = req.body;

  if (!smiles || !Array.isArray(smiles)) {
    return res.status(400).json({ error: "smiles array is required" });
  }

  // Ensure SDF directory exists
  if (!fs.existsSync(SDF_DIR)) {
    fs.mkdirSync(SDF_DIR, { recursive: true });
  }

  const sdfPaths = [];
  const errors = [];

  // Filter and process each chemical
  const validChemicals = smiles.filter(s => s && typeof s === 'string' && s.trim());

  const sdfPromises = validChemicals.map(s => {
    if (!s) return Promise.resolve();

    // Get processor info (SMILES first, crystal as backup)
    const { type, processor, backup } = getChemicalProcessor(s);
    
    // Check if SDF already exists (with smart filename detection)
    if (!overwrite) {
      const existingSdfPath = findGeneratedSdfFile(s, SDF_DIR);
      if (existingSdfPath) {
        console.log(`✅ Using existing file: ${s} → ${existingSdfPath}`);
        sdfPaths.push(existingSdfPath);
        return Promise.resolve();
      }
    }

    return new Promise((resolve, reject) => {
      console.log(`🧬 Generating structure for: ${s} (trying SMILES first)`);
      
      // Try SMILES first
      const { command, args } = processor(s, overwrite);
      const pythonProcess = spawn(command, args);
      
      let pythonOutput = '';
      pythonProcess.stdout.on('data', data => {
        const output = data.toString().trim();
        console.log(`Python Output: ${output}`);
        pythonOutput += output + '\n';
      });
      
      pythonProcess.stderr.on('data', data =>
        console.error(`Error: ${data.toString().trim()}`));

      pythonProcess.on('close', code => {
        if (code === 0) {
          // Extract actual filename from Python output or check for existing files
          const actualSdfPath = findGeneratedSdfFile(s, SDF_DIR);
          if (actualSdfPath) {
            console.log(`✅ Successfully generated SMILES structure: ${s} → ${actualSdfPath}`);
            sdfPaths.push(actualSdfPath);
            resolve();
          } else {
            console.log(`⚠️ SMILES succeeded but couldn't find SDF file for ${s}, trying crystal backup...`);
            tryBackupCrystal();
          }
        } else {
          // SMILES failed, try crystal/mineral as backup
          console.log(`⚠️ SMILES failed for ${s}, trying crystal/mineral backup...`);
          tryBackupCrystal();
        }
      });
      
      function tryBackupCrystal() {
        const { command: backupCmd, args: backupArgs } = backup(s, overwrite);
        const backupProcess = spawn(backupCmd, backupArgs);

        backupProcess.stdout.on('data', data =>
          console.log(`Backup Output: ${data.toString().trim()}`));
        backupProcess.stderr.on('data', data =>
          console.error(`Backup Error: ${data.toString().trim()}`));

        backupProcess.on('close', backupCode => {
          if (backupCode === 0) {
            const actualSdfPath = findGeneratedSdfFile(s, SDF_DIR);
            if (actualSdfPath) {
              console.log(`✅ Successfully generated crystal structure: ${s} → ${actualSdfPath}`);
              sdfPaths.push(actualSdfPath);
              resolve();
            } else {
              const errorMsg = `Crystal generation succeeded but couldn't find SDF file for ${s}`;
              console.error(errorMsg);
              errors.push(errorMsg);
              reject(new Error(errorMsg));
            }
          } else {
            const errorMsg = `Both SMILES and crystal generation failed for ${s}`;
            console.error(errorMsg);
            errors.push(errorMsg);
            reject(new Error(errorMsg));
          }
        });
      }
    });
  });

  await Promise.allSettled(sdfPromises);

  const response = { 
    sdfPaths: sdfPaths.filter(p => p),
    message: `Generated ${sdfPaths.length} 3D structures from ${validChemicals.length} chemicals`
  };

  if (errors.length > 0) {
    response.errors = errors;
    response.message += ` (${errors.length} failed)`;
  }

  res.json(response);
});

// ==================== SERVER STARTUP ====================
// Auto-deployment test - timestamp: 2025-07-09
// Only start servers in local development (NOT in Cloud Functions, Netlify, or tests)
const isCloudFunction = process.env.FUNCTION_NAME || process.env.FUNCTION_TARGET || process.env.K_SERVICE || process.env.GOOGLE_CLOUD_PROJECT || process.env.GCP_PROJECT;
const isNetlify = process.env.NETLIFY;
const isTestMode = process.env.NODE_ENV === 'test' || process.env.JEST_WORKER_ID;
const isServerless = isCloudFunction || isNetlify;

if (!isServerless && !isTestMode) {
  // Local development mode
  const localIP = getLocalIPAddress();

  app.listen(PORT, '0.0.0.0', () => {
    console.log(`🚀 Local server running on http://0.0.0.0:${PORT}`);
  });

  // Start HTTPS server only in local development
  if (process.env.NODE_ENV !== 'production') {
    const credentials = generateSelfSignedCert(localIP);
    if (credentials) {
      https.createServer(credentials, app).listen(HTTPS_PORT, '0.0.0.0', () => {
        console.log(`🔒 HTTPS server running on https://${localIP}:${HTTPS_PORT}`);
        console.log(`📱 Visit this on your phone browser: https://${localIP}:${HTTPS_PORT}`);
      });
    } else {
      console.log("⚠️ HTTPS not started — use ngrok if needed");
    }
  }
} else {
  // Serverless mode - server handled by platform
  if (isCloudFunction) {
    console.log(`☁️ Running in Cloud Functions mode - server handled by platform`);
  } else if (isNetlify) {
    console.log(`◈ Running in Netlify mode - server handled by platform`);
  }
  
  // Check for required environment variables in production
  if (process.env.NODE_ENV === 'production' && !process.env.OPENAI_API_KEY) {
    console.error('❌ OPENAI_API_KEY environment variable is required for production');
    process.exit(1);
  }
}

// Always export the app for Cloud Functions
module.exports = app;