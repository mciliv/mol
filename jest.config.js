module.exports = {
  testEnvironment: 'node',
  testMatch: [
    '**/test/smoke.test.js'
  ],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**',
    '!**/coverage/**'
  ],
  testTimeout: 10000,
  verbose: true
}; 