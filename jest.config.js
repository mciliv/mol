module.exports = {
  testEnvironment: 'node',
  setupFilesAfterEnv: ['<rootDir>/test/fixtures/setup.js'],
  testMatch: [
    '<rootDir>/test/**/*.test.js',
    '<rootDir>/test/**/test_*.py'
  ],
  moduleDirectories: ['node_modules'],
  testPathIgnorePatterns: ['/node_modules/'],
  collectCoverageFrom: [
    'backend/**/*.js',
    'frontend/**/*.js',
    '!**/node_modules/**'
  ]
}; 