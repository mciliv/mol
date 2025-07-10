const OpenAI = require("openai");
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
- "chemicals": dictionary with "smiles" key containing an array of constituent molecules

IMPORTANT: 
- Include any constituents as SMILES array (use proper SMILES notation, not chemical formulas)
- Break down materials into their molecular components when possible
- Use SMILES format like "O" for water, "CCO" for ethanol, NOT formulas like "H2O" or "C2H6O"
- If you see a person or body part, return generic human body composition (water, proteins, etc.)

Chemical types and their representations:
${instructions}

Example responses:
{
  "object": "A glass of water",
  "chemicals": {
    "smiles": ["O"]
  }
}

{
  "object": "Limestone (calcium carbonate)",
  "chemicals": {
    "smiles": ["C(=O)([O-])[O-].[Ca+2]"]
  }
}

{
  "object": "Alcoholic beverage",
  "chemicals": {
    "smiles": ["CCO", "O"]
  }
}

{
  "object": "Plastic bottle (PET)",
  "chemicals": {
    "smiles": ["O=C(C1=CC=C(CO)C=C1)OC"]
  }
}

{
  "object": "Human hand",
  "chemicals": {
    "smiles": ["O", "NCCNCCNCCN", "CC(C)CC(N)C(=O)O"]
  }
}

{
  "object": "Coffee cup (ceramic)",
  "chemicals": {
    "smiles": ["O"]
  }
}`;
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
        chemicals: parsed.chemicals || {}
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
        chemicals: parsed.chemicals || {}
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
          chemicals: {
            smiles: ["O", "NCCNCCNCCN", "CC(C)CC(N)C(=O)O", "CCCCCCCCCCCCCCCC(=O)O"]
          }
        };
      }
      
      return {
        object: "Analysis failed",
        chemicals: {
          smiles: []
        }
      };
    }
  }
}

module.exports = AIAnalyzer; 