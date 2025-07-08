process.env.NODE_ENV ||= 'development';

if (process.env.NODE_ENV !== 'production') {
  try {
    require('dotenv').config();
  } catch (err) {
    console.log('⚠️ dotenv not found, using environment variables directly');
  }
}

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

const app = express();
const PORT = process.env.PORT || 8080;
const HTTPS_PORT = process.env.HTTPS_PORT || 3001;
const SDF_DIR = path.join(__dirname, "sdf_files");

if (!fs.existsSync(SDF_DIR)) fs.mkdirSync(SDF_DIR, { recursive: true });

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

function sdf(s, overwrite) {
  let command = "python";
  let args = ["sdf.py", s, "--dir", SDF_DIR];
  if (overwrite) args.push("--overwrite");
  return { command, args };
}

async function getSmilesForObject(object, imageBase64 = null, croppedImageBase64 = null) {
  const client = new OpenAI({ apiKey: process.env["OPENAI_API_KEY"] });
  
  let text = `Identify: ${object}. List as many relevant chemical structures as SMILES strings. Consider the specific material, compound, or substance shown.`;
  
  let content = [{ text: text, type: "input_text" }];
  
  if (imageBase64) {
    console.log(`🖼️  Using image context for enhanced SMILES analysis`);
    content.push({ detail: "high", type: "input_image", image_url: `data:image/jpeg;base64,${imageBase64}` });
    text = `Analyze the object "${object}" in this image. List relevant chemical structures as SMILES strings. Use visual details to be more specific about materials, compounds, or substances shown.`;
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
  
  if (!parsed.output_parsed || !parsed.output_parsed.smiles) {
    console.log("⚠️ No SMILES found in OpenAI response");
    return [];
  }
  
  return parsed.output_parsed.smiles;
}

async function identifyObjectFromImage(imageBase64, croppedImageBase64, x, y) {
  const client = new OpenAI({ apiKey: process.env["OPENAI_API_KEY"] });
  
  const identificationText = `You are a chemical analysis expert. The user clicked at coordinate (X: ${Math.round(x)}, Y: ${Math.round(y)}) in this image. Identify the specific substance, material, or chemical compound at that location. Be as chemically specific as possible - consider molecular composition, material type, active ingredients, or chemical structure when naming the object.`;

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

if (process.env.NODE_ENV === "development") {
  const livereload = require("livereload");
  const connectLivereload = require("connect-livereload");

  try {
    const liveReloadServer = livereload.createServer();
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

app.use((req, res, next) => {
  console.log(`Incoming request: ${req.method} ${req.url}`);
  next();
});

app.use(cors());
app.use(express.json({ limit: "50mb" }));
app.use(express.static(__dirname));
app.use("/sdf_files", express.static(SDF_DIR));
app.use("/favicon.ico", express.static(path.join(__dirname, "favicon.ico")));

app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "index.html"));
});

app.post("/image-molecules", async (req, res) => {
  const { imageBase64, croppedImageBase64, x, y } = req.body;

  if (!imageBase64) {
    return res.status(400).json({ error: "No image data provided" });
  }

  try {
    const identifiedObject = await identifyObjectFromImage(imageBase64, croppedImageBase64, x, y);
    console.log(`🔍 Identified object: ${identifiedObject}`);

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

app.post('/generate-sdfs', async (req, res) => {
  const { smiles, overwrite = false } = req.body;

  if (!smiles || !Array.isArray(smiles)) {
    return res.status(400).json({ error: "smiles array is required" });
  }

  if (!fs.existsSync(SDF_DIR)) {
    fs.mkdirSync(SDF_DIR, { recursive: true });
  }

  const sdfPaths = [];
  const errors = [];

  const validSmiles = smiles.filter(s => s && typeof s === 'string' && s.trim());

  const sdfPromises = validSmiles.map(s => {
    if (!s) return Promise.resolve();

    const sdfPath = `/sdf_files/${s}.sdf`;
    const fullPath = path.join(SDF_DIR, `${s}.sdf`);

    if (fs.existsSync(fullPath) && !overwrite) {
      sdfPaths.push(sdfPath);
      return Promise.resolve();
    }

    return new Promise((resolve, reject) => {
      const { command, args } = sdf(s, overwrite);
      const pythonProcess = spawn(command, args);

      pythonProcess.stdout.on('data', data =>
        console.log(`Python Output: ${data.toString().trim()}`));
      pythonProcess.stderr.on('data', data =>
        console.error(`Error: ${data.toString().trim()}`));

      pythonProcess.on('close', code => {
        if (code === 0) {
          sdfPaths.push(sdfPath);
          resolve();
        } else {
          const errorMsg = `SDF generation failed for ${s}`;
          console.error(errorMsg);
          errors.push(errorMsg);
          reject(new Error(errorMsg));
        }
      });
    });
  });

  await Promise.allSettled(sdfPromises);

  if (errors.length > 0) {
    console.error(`SDF generation errors: ${errors.join(', ')}`);
    return res.json({ 
      sdfPaths: sdfPaths.filter(p => p), 
      errors,
      message: `Generated ${sdfPaths.length} of ${validSmiles.length} SDFs` 
    });
  }

  res.json({ message: "Files generated", sdfPaths });
});

const isCloudFunction = process.env.FUNCTION_NAME || process.env.FUNCTION_TARGET || process.env.K_SERVICE || process.env.GOOGLE_CLOUD_PROJECT;
const isVercel = process.env.VERCEL || process.env.VERCEL_ENV;
const isNetlify = process.env.NETLIFY;
const isRailway = process.env.RAILWAY_ENVIRONMENT;
const isTestMode = process.env.NODE_ENV === 'test' || process.env.JEST_WORKER_ID;
const isServerless = isCloudFunction || isVercel || isNetlify || isRailway;

if (!isServerless && !isTestMode) {
  const localIP = getLocalIPAddress();

  app.listen(PORT, '0.0.0.0', () => {
    console.log(`🚀 Local server running on http://0.0.0.0:${PORT}`);
  });

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
  if (isCloudFunction) {
    console.log(`☁️ Running in Cloud Functions mode - server handled by platform`);
  } else if (isVercel) {
    console.log(`▲ Running in Vercel mode - server handled by platform`);
  } else if (isNetlify) {
    console.log(`◈ Running in Netlify mode - server handled by platform`);
  } else if (isRailway) {
    console.log(`🚂 Running in Railway mode - server handled by platform`);
  }
  
  if (process.env.NODE_ENV === 'production' && !process.env.OPENAI_API_KEY) {
    console.error('❌ OPENAI_API_KEY environment variable is required for production');
    process.exit(1);
  }
}

module.exports = app;