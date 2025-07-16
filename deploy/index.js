const functions = require('@google-cloud/functions-framework');
const app = require('./server.js');

// Register the Express app as a Cloud Function
functions.http('molecularAnalysis', app);

console.log('🟢 Google Cloud Functions - Molecular Analysis App Ready'); 