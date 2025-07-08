const { z } = require("zod");

// Schema for object identification only (from images)
const ObjectIdentificationSchema = z.object({
  object: z.string(),
});

// Schema for image-based molecule analysis (includes object identification)
const ImageMoleculeSchema = z.object({
  object: z.string(),
  smiles: z.array(z.string()),
});

// Schema for text-based molecule analysis (just SMILES list)
const TextMoleculeSchema = z.object({
  smiles: z.array(z.string()),
});

// Schema for request body in /list-molecules-text
const ListMoleculesTextRequestSchema = z.object({
  object: z.string(),
});

module.exports = {
  ObjectIdentificationSchema,
  ImageMoleculeSchema,
  TextMoleculeSchema,
  ListMoleculesTextRequestSchema,
};
