const functions = require('@google-cloud/functions-framework');
const app = require('./server.js');

// Register HTTP function for Google Cloud Functions
functions.http('molecularAnalysis', (req, res) => {
  // Handle the request using the Express app
  app(req, res);
});

console.log('🟢 Google Cloud Functions - Molecular Analysis App Ready'); 