// Server-side code only
const express = require("express");
const cors = require("cors");
const fs = require("fs");
const path = require("path");
const { spawn } = require("child_process");
const { env } = require("process");
const OpenAI = require("openai");
const {
  ImageMoleculeSchema,
  TextMoleculeSchema,
  ListMoleculesTextRequestSchema,
} = require("./schemas");


const app = express();
const PORT = 3000;
const SDF_DIR = path.join(__dirname, "sdf_files");


// Ensure SDF directory exists
if (!fs.existsSync(SDF_DIR)) {
  fs.mkdirSync(SDF_DIR, { recursive: true });
}


// Middleware
app.use((req, res, next) => {
  console.log(`Incoming request: ${req.method} ${req.url}`);
  next();
});


app.use(cors());
app.use(express.json({ limit: "50mb" }));
app.use(express.static(__dirname));
app.use("/sdf_files", express.static(SDF_DIR));
app.use("/favicon.ico", express.static(path.join(__dirname, "favicon.ico")));


// Routes
app.get("/", (req, res) => {
  res.sendFile(path.join(__dirname, "index.html"));
});


app.post("/list-molecules", async (req, res) => {
  const { imageBase64, croppedImageBase64, x, y } = req.body;


  if (!imageBase64) {
    return res.status(400).json({ error: "No image data provided" });
  }


  try {
    const client = new OpenAI({
      apiKey: process.env["OPENAI_API_KEY"],
    });


    const text = `The goal of this application is to identify a substance/material, and the user clues us as to what they're identifying by clicking on the object in the image. So, what is at coordinate (X: ${Math.round(x)}, Y: ${Math.round(y)}) in the image?`;


    fs.writeFileSync("image.jpg", imageBase64, "base64");
    fs.writeFileSync("cropped_image.jpg", croppedImageBase64, "base64");


    const response = await client.responses.create({
      input: [
        {
          content: [
            {
              text: text,
              type: "input_text",
            },
            {
              detail: "high",


              type: "input_image",
              image_url: `data:image/jpeg;base64,${imageBase64}`,
            },
            {
              detail: "high",
              type: "input_image",
              image_url: `data:image/jpeg;base64,${croppedImageBase64}`,
            },
          ],
          role: "user",
        },
      ],
      model: "gpt-4.1",
      text: {
        format: {
          type: "json_schema",
          schema: ImageMoleculeSchema.shape,
        },
      },
    });


    const validatedResponse = ImageMoleculeSchema.parse(
      response.results[0].output,
    );
    return res.json({ output: validatedResponse });
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


// Route for manual object input to list molecule SMILES
app.post("/list-molecules-text", async (req, res) => {
  try {
    const validatedData = ListMoleculesTextRequestSchema.parse(req.body);
    const { object } = validatedData;


    const client = new OpenAI({ apiKey: process.env["OPENAI_API_KEY"] });
    const text = `The user entered: ${object}. List as many relevant chemical structures as SMILES in a direct JSON array without any markdown or additional text.`;


    const aiRes = await client.responses.create({
      input: [{ content: [{ text: text, type: "input_text" }], role: "user" }],
      model: "gpt-4.1",
      text: {
        format: { type: "json_schema", schema: TextMoleculeSchema.shape },
      },
    });


    const validatedResponse = TextMoleculeSchema.parse(aiRes.results[0].output);
    return res.json({ output: validatedResponse });
  } catch (error) {
    console.error("Error in /list-molecules-text handler:", error);
    if (error.name === "ZodError") {
      return res
        .status(400)
        .json({ error: "Invalid request data", details: error.errors });
    }
    res.status(500).json({ error: error.message });
  }
});


// Helper functions
function sdf(s, overwrite) {
  let command = "python";
  let args = ["sdf.py", s, "--dir", SDF_DIR];
  if (overwrite) args.push("--overwrite");
  return { command, args };
}


// Start server
app.listen(PORT, () =>
  console.log(`Node server running on http://localhost:${PORT}`),
);


module.exports = app;
