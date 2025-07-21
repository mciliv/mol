// test/fixtures/setup.js - Jest setup file for mocking dependencies

// Mock OpenAI before any modules are imported
jest.mock("openai", () => {
  return {
    OpenAI: jest.fn().mockImplementation(() => ({
      responses: {
        parse: jest.fn().mockResolvedValue({
          output_parsed: {
            object: "test object",
            smiles: ["O", "CCO"],
          },
        }),
      },
    })),
  };
});

// Set test environment variables
process.env.NODE_ENV = "test";
process.env.OPENAI_API_KEY = "test-key";
