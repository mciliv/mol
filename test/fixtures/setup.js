// Test environment setup
const { TextEncoder, TextDecoder } = require('util');

// Add TextEncoder/TextDecoder globally for tests
global.TextEncoder = TextEncoder;
global.TextDecoder = TextDecoder;

// Mock console methods to reduce noise in tests
const originalConsoleLog = console.log;
const originalConsoleError = console.error;

// Only show important test output
console.log = (...args) => {
  const message = args.join(' ');
  if (message.includes('✅') || message.includes('❌') || message.includes('🧪')) {
    originalConsoleLog(...args);
  }
};

console.error = (...args) => {
  const message = args.join(' ');
  if (!message.includes('Database connection failed') && !message.includes('PostgreSQL')) {
    originalConsoleError(...args);
  }
};

// Increase timeout for integration tests
jest.setTimeout(30000);

// Mock WebSocket if running in test environment
if (typeof WebSocket === 'undefined') {
  global.WebSocket = class MockWebSocket {
    constructor() {
      this.readyState = 1; // OPEN
    }
    send() {}
    close() {}
  };
}
