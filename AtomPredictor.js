const { OpenAI } = require("openai");
const { ObjectIdentificationSchema, CHEMICAL_REPRESENTATIONS } = require("./schemas");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ apiKey });
    this.chemicalInstructions = this.buildChemicalInstructions();
  }

  buildChemicalInstructions() {
    const instructions = Object.entries(CHEMICAL_REPRESENTATIONS)
      .map(([type, description]) => `${type}: ${description}`)
      .join('\n');
    
    return `Analyze the image and identify what the user clicked on. Return a JSON object with:
- "object": brief description of what you see  
- "chemicals": array of objects with "name" and "smiles" fields (include both major and minor components when possible)

IMPORTANT:
- Focus on the EXACT pixel location the user clicked on, not the general area or major objects in the image
- The red box in the cropped image marks the precise pixel coordinates - analyze what's at that specific point
- Be as chemically specific and realistic as possible. For food, beverages, and natural materials, include both major and minor chemical constituents (e.g., minerals, ions, organic acids, polyphenols, sugars, etc.).
- CRITICAL: Use proper SMILES notation ONLY, never molecular formulas like C3N2O2H6 or trivial names.
- SMILES examples: "CCO" for ethanol, "C(C(=O)O)N" for glycine, "C1=CC=CC=C1" for benzene
- NEVER use molecular formulas like "C2H6O" or "C3N2O2H6" - these are not SMILES
- WRONG: {"name": "Histidine", "smiles": "C3N2O2H6"} - this is a molecular formula
- CORRECT: {"name": "Histidine", "smiles": "C1=C(NC=N1)CC(C(=O)O)N"} - this is proper SMILES
- For mixtures (e.g., wine, seawater), list the most abundant and characteristic molecules.
- For generic objects (e.g., "alcoholic beverage"), prefer to infer the likely type (e.g., wine, beer, spirits) and break down accordingly.
- If you see a person or body part, return a realistic set of major body constituents (water, proteins, lipids, etc.).
- For complex biological structures (collagen, DNA, proteins, polymers), provide the actual complex molecular structures when possible, not just simple components.
- For materials like plastics, fibers, or biological tissues, include the characteristic polymer or macromolecular structures.
- When you see biological tissues like skin, hair, nails, tendons, or connective tissue, these contain COLLAGEN and KERATIN - return the actual protein components.
- COLLAGEN contains: glycine, proline, hydroxyproline, and collagen peptide segments
- KERATIN contains: cysteine, serine, glycine, and keratin protein segments
- Do NOT return random small molecules when analyzing biological proteins - identify the actual protein components.

Chemical types and their representations:
${instructions}

Example responses:
{
  "object": "Tap water",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Sodium ion", "smiles": "[Na+]"},
    {"name": "Chloride ion", "smiles": "[Cl-]"},
    {"name": "Calcium ion", "smiles": "[Ca+2]"},
    {"name": "Magnesium ion", "smiles": "[Mg+2]"},
    {"name": "Potassium ion", "smiles": "[K+]"},
    {"name": "Bicarbonate ion", "smiles": "[HCO3-]"}
  ]
}

{
  "object": "Mineral water",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Sodium ion", "smiles": "[Na+]"},
    {"name": "Chloride ion", "smiles": "[Cl-]"},
    {"name": "Calcium ion", "smiles": "[Ca+2]"},
    {"name": "Magnesium ion", "smiles": "[Mg+2]"},
    {"name": "Potassium ion", "smiles": "[K+]"},
    {"name": "Bicarbonate ion", "smiles": "[HCO3-]"},
    {"name": "Sulfate ion", "smiles": "[SO4-2]"},
    {"name": "Iron(II) ion", "smiles": "[Fe+2]"},
    {"name": "Zinc ion", "smiles": "[Zn+2]"}
  ]
}

{
  "object": "Seawater",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Sodium ion", "smiles": "[Na+]"},
    {"name": "Chloride ion", "smiles": "[Cl-]"},
    {"name": "Magnesium ion", "smiles": "[Mg+2]"},
    {"name": "Sulfate ion", "smiles": "[SO4-2]"},
    {"name": "Calcium ion", "smiles": "[Ca+2]"},
    {"name": "Potassium ion", "smiles": "[K+]"},
    {"name": "Bicarbonate ion", "smiles": "[HCO3-]"},
    {"name": "Bromide ion", "smiles": "[Br-]"},
    {"name": "Strontium ion", "smiles": "[Sr+2]"},
    {"name": "Fluoride ion", "smiles": "[F-]"}
  ]
}

{
  "object": "Wine",
  "chemicals": [
    {"name": "Ethanol", "smiles": "CCO"},
    {"name": "Water", "smiles": "O"},
    {"name": "Tartaric acid", "smiles": "OC(C(O)C(O)=O)C(O)=O"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Fructose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Malic acid", "smiles": "C(C(=O)O)C(C(=O)O)O"},
    {"name": "Citric acid", "smiles": "C(C(=O)O)C(O)(C(=O)O)O"},
    {"name": "Resveratrol", "smiles": "C1=CC(=C(C=C1)O)O"},
    {"name": "Phenol", "smiles": "C1=CC=C(C=C1)O"},
    {"name": "Benzyl alcohol", "smiles": "C1=CC=C(C=C1)CO"},
    {"name": "Benzoic acid", "smiles": "C1=CC=C(C=C1)C(=O)O"}
  ]
}

{
  "object": "Beer",
  "chemicals": [
    {"name": "Ethanol", "smiles": "CCO"},
    {"name": "Water", "smiles": "O"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Fructose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Acetic acid", "smiles": "CC(=O)O"},
    {"name": "Pyruvic acid", "smiles": "CC(=O)C(=O)O"},
    {"name": "Phenol", "smiles": "C1=CC=C(C=C1)O"},
    {"name": "Benzyl alcohol", "smiles": "C1=CC=C(C=C1)CO"},
    {"name": "Benzoic acid", "smiles": "C1=CC=C(C=C1)C(=O)O"}
  ]
}

{
  "object": "Coffee",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Caffeine", "smiles": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"},
    {"name": "Chlorogenic acid", "smiles": "C1=CC(=C(C=C1)O)O"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Phenol", "smiles": "C1=CC=C(C=C1)O"},
    {"name": "Benzyl alcohol", "smiles": "C1=CC=C(C=C1)CO"}
  ]
}

{
  "object": "Plastic bottle (PET)",
  "chemicals": [
    {"name": "PET repeat unit", "smiles": "O=C(C1=CC=CC=C1)OC2=CC=CC=C2C(=O)O"},
    {"name": "Isobutyric acid", "smiles": "CC(C)C(=O)O"},
    {"name": "Acetic acid", "smiles": "CC(=O)O"},
    {"name": "Ethanol", "smiles": "CCO"}
  ]
}

{
  "object": "Human hand",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Leucine", "smiles": "CC(C)C(C(=O)O)N"},
    {"name": "Glycine", "smiles": "C(C(=O)O)N"},
    {"name": "Palmitic acid", "smiles": "CCCCCCCCCCCCCCCC(=O)O"},
    {"name": "Lactic acid", "smiles": "C(C(=O)O)O"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"}
  ]
}

{
  "object": "Grape",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Fructose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"},
    {"name": "Tartaric acid", "smiles": "OC(C(O)C(O)=O)C(O)=O"},
    {"name": "Malic acid", "smiles": "C(C(=O)O)C(C(=O)O)O"},
    {"name": "Resveratrol", "smiles": "C1=CC(=C(C=C1)O)O"}
  ]
}

{
  "object": "Collagen fiber",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Glycine", "smiles": "C(C(=O)O)N"},
    {"name": "Proline", "smiles": "C1CC(NC1)C(=O)O"},
    {"name": "Hydroxyproline", "smiles": "C1C(C(NC1)C(=O)O)O"},
    {"name": "Alanine", "smiles": "CC(C(=O)O)N"},
    {"name": "Collagen peptide segment", "smiles": "NC(C(=O)NC(C(=O)NC1CCC(N1)C(=O)O)C)C"}
  ]
}

{
  "object": "Plastic material",
  "chemicals": [
    {"name": "Polyethylene segment", "smiles": "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC"},
    {"name": "Polypropylene segment", "smiles": "CC(C)CC(C)CC(C)CC(C)CC(C)CC(C)C"},
    {"name": "PET polymer segment", "smiles": "O=C(C1=CC=C(C=C1)C(=O)OCCOC(=O)C2=CC=C(C=C2)C(=O)OCCOC(=O)C3=CC=C(C=C3)C(=O)O)O"}
  ]
}

{
  "object": "Hair strand (keratin)",
  "chemicals": [
    {"name": "Water", "smiles": "O"},
    {"name": "Cysteine", "smiles": "C(C(C(=O)O)N)S"},
    {"name": "Serine", "smiles": "C(C(C(=O)O)N)O"},
    {"name": "Glycine", "smiles": "C(C(=O)O)N"},
    {"name": "Alanine", "smiles": "CC(C(=O)O)N"},
    {"name": "Keratin peptide segment", "smiles": "NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)O)CS)CO)C)CS"}
  ]
}
`;
  }

  async analyzeImage(imageBase64, croppedImageBase64 = null, x = null, y = null, cropMiddleX = null, cropMiddleY = null, cropSize = null) {
    try {
      const messages = [
        {
          role: "user",
          content: [
            {
              type: "text",
              text: this.chemicalInstructions
            },
            {
              type: "image_url",
              image_url: {
                url: `data:image/jpeg;base64,${imageBase64}`,
                detail: "high"
              }
            }
          ]
        }
      ];

      // Add cropped region if available
      if (croppedImageBase64) {
        let focusText = `CRITICAL: Focus on the EXACT pixel at coordinates (${x}, ${y}) in the original image. The user clicked on this specific point, not the general area. Here's a cropped view:`;
        
        // Add information about the middle pixel if available
        if (cropMiddleX !== null && cropMiddleY !== null && cropSize !== null) {
          focusText += ` The red box in the cropped image marks the exact middle pixel (${cropMiddleX}, ${cropMiddleY}) of the ${cropSize}x${cropSize} crop region. This is the precise point the user clicked on. Analyze ONLY what's at this specific pixel location.`;
        }
        
        messages[0].content.push({
          type: "text",
          text: focusText
        });
        messages[0].content.push({
          type: "image_url",
          image_url: {
            url: `data:image/jpeg;base64,${croppedImageBase64}`,
            detail: "high"
          }
        });
      }

      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages,
        max_tokens: 1000,
        temperature: 0.1
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);
      
      // Extract SMILES for backward compatibility
      const smiles = parsed.chemicals ? parsed.chemicals.map(chem => chem.smiles).filter(Boolean) : [];
      
      return {
        object: parsed.object || "Unknown object",
        chemicals: parsed.chemicals || [],
        smiles: smiles  // Add for backward compatibility with tests
      };

    } catch (error) {
      console.error("AI analysis error:", error);
      throw new Error(`AI analysis failed: ${error.message}`);
    }
  }

  async analyzeText(object) {
    try {
      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages: [
          {
            role: "user",
            content: `Analyze this object: "${object}". ${this.chemicalInstructions}`
          }
        ],
        max_tokens: 1000,
        temperature: 0.1
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);
      
      // Extract SMILES for backward compatibility
      const smiles = parsed.chemicals ? parsed.chemicals.map(chem => chem.smiles).filter(Boolean) : [];
      
      return {
        object: parsed.object || object,
        chemicals: parsed.chemicals || [],
        smiles: smiles  // Add for backward compatibility with tests
      };

    } catch (error) {
      console.error("AI text analysis error:", error);
      throw new Error(`AI text analysis failed: ${error.message}`);
    }
  }

  parseAIResponse(content) {
    try {
      // Try to extract JSON from the response
      const jsonMatch = content.match(/\{[\s\S]*\}/);
      if (jsonMatch) {
        return JSON.parse(jsonMatch[0]);
      }
      
      // Fallback: try to parse the entire content as JSON
      return JSON.parse(content);
    } catch (error) {
      console.error("Failed to parse AI response:", content);
      
      // Check if AI refused to analyze people - provide generic human body composition
      if (content.includes("unable to identify or analyze people")) {
        return {
          object: "Human body (generic composition)",
          chemicals: [
            {"name": "Water", "smiles": "O"},
            {"name": "Polyethylene glycol", "smiles": "NCCNCCNCCN"},
            {"name": "Leucine", "smiles": "CC(C)CC(N)C(=O)O"},
            {"name": "Palmitic acid", "smiles": "CCCCCCCCCCCCCCCC(=O)O"}
          ]
        };
      }
      
      return {
        object: "Analysis failed",
        chemicals: []
      };
    }
  }
}

module.exports = AtomPredictor; 