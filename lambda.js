const serverlessExpress = require('@vendia/serverless-express');
const app = require('./server.js');

// Create serverless Express instance
const serverlessExpressInstance = serverlessExpress({ app });

// Lambda handler function
module.exports.handler = async (event, context) => {
  console.log('🚀 AWS Lambda - Molecular Analysis Request');
  console.log('Event:', JSON.stringify(event, null, 2));
  
  // Handle the request through serverless Express
  const result = await serverlessExpressInstance(event, context);
  
  console.log('✅ Response sent');
  return result;
}; 