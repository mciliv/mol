const { OpenAI } = require("openai");
const { ObjectIdentificationSchema, CHEMICAL_REPRESENTATIONS } = require("./schemas");

class AIAnalyzer {
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
- Be as chemically specific and realistic as possible. For food, beverages, and natural materials, include both major and minor chemical constituents (e.g., minerals, ions, organic acids, polyphenols, sugars, etc.).
- Use proper SMILES notation, not chemical formulas or trivial names.
- For mixtures (e.g., wine, seawater), list the most abundant and characteristic molecules.
- For generic objects (e.g., "alcoholic beverage"), prefer to infer the likely type (e.g., wine, beer, spirits) and break down accordingly.
- If you see a person or body part, return a realistic set of major body constituents (water, proteins, lipids, etc.).

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
`;
  }

  async analyzeImage(imageBase64, croppedImageBase64 = null, x = null, y = null) {
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
        messages[0].content.push({
          type: "text",
          text: `Focus on the region around coordinates (${x}, ${y}). Here's a cropped view:`
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
      
      return {
        object: parsed.object || "Unknown object",
        chemicals: parsed.chemicals || []
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
      
      return {
        object: parsed.object || object,
        chemicals: parsed.chemicals || []
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

module.exports = AIAnalyzer; 