const { z } = require("zod/v4");


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
  ImageMoleculeSchema,
  TextMoleculeSchema,
  ListMoleculesTextRequestSchema,
};
