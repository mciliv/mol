const { OpenAI } = require("openai");
const {
  ObjectIdentificationSchema,
  CHEMICAL_REPRESENTATIONS,
} = require("../schemas/schemas");

class AtomPredictor {
  constructor(apiKey) {
    this.client = new OpenAI({ 
      apiKey,
      timeout: 30000, // 30 seconds timeout
      maxRetries: 2   // Retry failed requests twice
    });
    this.chemicalInstructions = this.buildChemicalInstructions();
  }

  buildChemicalInstructions() {
    return `Analyze the object and provide a JSON response with chemical components.

Response format:
{
  "object": "Object name",
  "chemicals": [
    {"name": "Chemical name", "smiles": "SMILES notation"}
  ]
}

Use standard SMILES notation.`;
  }

  async analyzeImage(
    imageBase64,
    croppedImageBase64 = null,
    x = null,
    y = null,
    cropMiddleX = null,
    cropMiddleY = null,
    cropSize = null,
  ) {
    try {
      const messages = [
        {
          role: "user",
          content: [
            {
              type: "text",
              text: this.chemicalInstructions,
            },
            {
              type: "image_url",
              image_url: {
                url: `data:image/jpeg;base64,${imageBase64}`,
                detail: "high",
              },
            },
          ],
        },
      ];

      // Add cropped region if available
      if (croppedImageBase64) {
        let focusText = `Here's a cropped view of the area of interest. Analyze the chemical composition of the material or substance visible in this region:`;

        messages[0].content.push({
          type: "text",
          text: focusText,
        });
        messages[0].content.push({
          type: "image_url",
          image_url: {
            url: `data:image/jpeg;base64,${croppedImageBase64}`,
            detail: "high",
          },
        });
      }

      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages,
        max_tokens: 1000,
        temperature: 0.1,
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);

      return {
        object: parsed.object || "Unknown object",
        chemicals: parsed.chemicals || [],
      };
    } catch (error) {
      console.error("AI analysis error:", error);
      
      // Enhanced error handling with specific error types
      let errorMessage = `AI analysis failed: ${error.message}`;
      
      if (error.code === 'ENOTFOUND' || error.code === 'ECONNREFUSED') {
        errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      } else if (error.name === 'APIConnectionError') {
        errorMessage = 'Unable to connect to AI service. Please check your internet connection.';
      } else if (error.status === 401) {
        errorMessage = 'API authentication failed. Please check your API key configuration.';
      } else if (error.status === 429) {
        errorMessage = 'Rate limit exceeded. Please wait a moment and try again.';
      } else if (error.status === 503) {
        errorMessage = 'AI service temporarily unavailable. Please try again in a few moments.';
      }
      
      throw new Error(errorMessage);
    }
  }

  async analyzeText(object) {
    try {
      const response = await this.client.chat.completions.create({
        model: "gpt-4o",
        messages: [
          {
            role: "user",
            content: `Analyze this object: "${object}". ${this.chemicalInstructions}`,
          },
        ],
        max_tokens: 1000,
        temperature: 0.1,
      });

      const content = response.choices[0].message.content;
      const parsed = this.parseAIResponse(content);

      return {
        object: parsed.object || object,
        chemicals: parsed.chemicals || [],
      };
    } catch (error) {
      console.error("AI text analysis error:", error);
      
      // Enhanced error handling with specific error types
      let errorMessage = `AI text analysis failed: ${error.message}`;
      
      if (error.code === 'ENOTFOUND' || error.code === 'ECONNREFUSED') {
        errorMessage = 'Network connection failed. Please check your internet connection and try again.';
      } else if (error.name === 'APIConnectionError') {
        errorMessage = 'Unable to connect to AI service. Please check your internet connection.';
      } else if (error.status === 401) {
        errorMessage = 'API authentication failed. Please check your API key configuration.';
      } else if (error.status === 429) {
        errorMessage = 'Rate limit exceeded. Please wait a moment and try again.';
      } else if (error.status === 503) {
        errorMessage = 'AI service temporarily unavailable. Please try again in a few moments.';
      }
      
      throw new Error(errorMessage);
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
      
      // Simple fallback - just return empty result
      return {
        object: "Unknown object",
        chemicals: []
      };
    }
  }
}

module.exports = AtomPredictor;
