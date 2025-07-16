const { z } = require("zod");

// Chemical representations - focus on SMILES constituents
const CHEMICAL_REPRESENTATIONS = {
  smiles: "Array of SMILES notation for constituent molecules (e.g., ['CCO', 'O'] for alcoholic beverage, ['C1=CC=CC=C1'] for benzene)"
};

// Schema for object identification with chemical names and SMILES
const ObjectIdentificationSchema = z.object({
  object: z.string().describe("The identified object or material"),
  chemicals: z.array(z.object({
    name: z.string().describe("Chemical name or common name"),
    smiles: z.string().describe("SMILES notation for the molecule")
  })).describe("Array of chemical constituents with names and SMILES notation")
});

// Schema for image-based molecular analysis
const ImageMoleculeSchema = z.object({
  imageBase64: z.string().describe("Base64 encoded image data"),
  croppedImageBase64: z.string().optional().describe("Base64 encoded cropped region"),
  x: z.number().optional().describe("X coordinate of click"),
  y: z.number().optional().describe("Y coordinate of click")
});

// Schema for text-based molecular analysis
const TextMoleculeSchema = z.object({
  object: z.string().describe("Text description of object or material")
});

// Schema for generating SDF files
const SdfGenerationSchema = z.object({
  smiles: z.array(z.string()).describe("Array of SMILES notation to convert to SDF files"),
  overwrite: z.boolean().optional().default(false).describe("Whether to overwrite existing SDF files")
});

module.exports = {
  ObjectIdentificationSchema,
  ImageMoleculeSchema,
  TextMoleculeSchema,
  SdfGenerationSchema,
  CHEMICAL_REPRESENTATIONS
};
