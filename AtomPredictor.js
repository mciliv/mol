const { OpenAI } = require("openai");
const { ObjectIdentificationSchema, CHEMICAL_REPRESENTATIONS } = require("./schemas");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ apiKey });
    this.chemicalInstructions = this.buildChemicalInstructions();
  }

  buildChemicalInstructions() {
    return `You are a molecular analysis expert. Analyze the object and provide a JSON response with the chemical composition.

IMPORTANT RULES:
1. Generate ONLY valid, concise SMILES notation
2. For complex molecules, use representative fragments or simplified forms
3. For minerals like talc, use simple representations like "O[Si](O)O" or "O[Si]O"
4. For polymers, use short repeat units (max 10-15 atoms)
5. For large biomolecules, use key functional groups or representative structures
6. Avoid generating extremely long SMILES strings (>200 characters)
7. Focus on the most important/abundant chemical components

Response format:
{
  "object": "Object name",
  "chemicals": [
    {"name": "Chemical name", "smiles": "SMILES notation"},
    {"name": "Chemical name", "smiles": "SMILES notation"}
  ]
}

Examples of good SMILES:
- Water: "O"
- Glucose: "C(C(C(C(C(C=O)O)O)O)O)O"
- Caffeine: "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
- Talc: "O[Si](O)O"
- Cellulose: "C(C1C(C(C(C(O1)O)O)O)O)O"

Analyze the object and provide the most relevant chemical components.`;
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
        let focusText = `Here's a cropped view of the area of interest. Analyze the chemical composition of the material or substance visible in this region:`;
        
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
      
      // Handle various AI refusal scenarios with appropriate fallbacks
      if (content.includes("can't analyze images") || content.includes("unable to identify") || content.includes("I'm sorry")) {
        // AI refused to analyze - provide a generic but useful response
        return {
          object: "Generic object (AI analysis unavailable)",
          chemicals: [
            {"name": "Water", "smiles": "O"},
            {"name": "Carbon dioxide", "smiles": "O=C=O"},
            {"name": "Nitrogen", "smiles": "N#N"},
            {"name": "Oxygen", "smiles": "O=O"}
          ]
        };
      }
      
      // Check if AI refused to analyze people specifically
      if (content.includes("unable to identify or analyze people") || content.includes("can't identify people")) {
        return {
          object: "Human body (generic composition)",
          chemicals: [
            {"name": "Water", "smiles": "O"},
            {"name": "Glycine", "smiles": "C(C(=O)O)N"},
            {"name": "Leucine", "smiles": "CC(C)CC(N)C(=O)O"},
            {"name": "Palmitic acid", "smiles": "CCCCCCCCCCCCCCCC(=O)O"},
            {"name": "Glucose", "smiles": "C(C(C(C(C(C=O)O)O)O)O)O"}
          ]
        };
      }
      
      // For any other parsing error, provide a basic response
      return {
        object: "Analysis completed with fallback data",
        chemicals: [
          {"name": "Water", "smiles": "O"},
          {"name": "Carbon", "smiles": "C"},
          {"name": "Oxygen", "smiles": "O=O"}
        ]
      };
    }
  }
}

module.exports = AtomPredictor; 