module.exports = {
  testEnvironment: 'node',
  setupFilesAfterEnv: ['<rootDir>/testing/fixtures/setup.js'],
  testMatch: [
    '<rootDir>/testing/**/*.test.js',
    '<rootDir>/testing/**/test_*.py'
  ],
  moduleDirectories: ['node_modules'],
  testPathIgnorePatterns: ['/node_modules/'],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**'
  ]
}; 