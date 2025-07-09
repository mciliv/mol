#!/usr/bin/env node

const { spawn, execSync } = require('child_process');
const fs = require('fs');
const path = require('path');

console.log('🚀 Auto-deployment watcher started for Google Cloud Functions');
console.log('📁 Watching: server.js, schemas.js, *.py files, package.json');
console.log('☁️  Target: molecular-analysis function in us-central1');
console.log('⏰ Debounce: 5 seconds (multiple changes batched)');

let deployTimeout = null;
let isDeploying = false;

// Get OpenAI API key from environment
const apiKey = process.env.OPENAI_API_KEY;
if (!apiKey) {
  console.error('❌ OPENAI_API_KEY environment variable required');
  process.exit(1);
}

function deploy() {
  if (isDeploying) {
    console.log('⏳ Deployment already in progress, skipping...');
    return;
  }

  isDeploying = true;
  const timestamp = new Date().toLocaleTimeString();
  console.log(`\n🔄 [${timestamp}] Starting deployment...`);

  const deployCmd = [
    'gcloud', 'functions', 'deploy', 'molecular-analysis',
    '--gen2',
    '--runtime', 'nodejs20', 
    '--trigger-http',
    '--allow-unauthenticated',
    '--memory', '1GB',
    '--timeout', '540s',
    '--region', 'us-central1',
    '--entry-point', 'molecularAnalysis',
    '--source', '.',
    '--set-env-vars', `OPENAI_API_KEY=${apiKey}`,
    '--quiet'  // Suppress interactive prompts
  ];

  const deployProcess = spawn(deployCmd[0], deployCmd.slice(1), {
    stdio: ['pipe', 'pipe', 'pipe'],
    cwd: path.join(__dirname, '..')
  });

  let output = '';
  deployProcess.stdout.on('data', (data) => {
    output += data.toString();
  });

  deployProcess.stderr.on('data', (data) => {
    output += data.toString();
  });

  deployProcess.on('close', (code) => {
    const endTime = new Date().toLocaleTimeString();
    if (code === 0) {
      console.log(`✅ [${endTime}] Deployment successful!`);
      console.log('🌐 URL: https://us-central1-mol-analysis-app.cloudfunctions.net/molecular-analysis');
    } else {
      console.log(`❌ [${endTime}] Deployment failed (exit code: ${code})`);
      if (output.includes('revision')) {
        console.log('📋 Partial output:', output.slice(-200));
      }
    }
    isDeploying = false;
    console.log('👁️  Watching for changes...\n');
  });
}

function scheduleDeployment(filename) {
  const timestamp = new Date().toLocaleTimeString();
  console.log(`📝 [${timestamp}] Change detected: ${filename}`);
  
  if (deployTimeout) {
    clearTimeout(deployTimeout);
  }
  
  deployTimeout = setTimeout(() => {
    deploy();
  }, 5000); // 5 second debounce
}

// Watch key files
const filesToWatch = [
  'server.js',
  'schemas.js', 
  'index.js',
  'package.json',
  'sdf.py',
  'crystal.py'
];

filesToWatch.forEach(file => {
  const fullPath = path.join(__dirname, '..', file);
  if (fs.existsSync(fullPath)) {
    fs.watchFile(fullPath, { interval: 1000 }, (curr, prev) => {
      if (curr.mtime !== prev.mtime) {
        scheduleDeployment(file);
      }
    });
    console.log(`👁️  Watching: ${file}`);
  }
});

console.log('\n✅ Auto-deployment watcher ready!');
console.log('💡 Make changes to watched files to trigger deployment');
console.log('🛑 Press Ctrl+C to stop\n');

// Keep the process alive
process.on('SIGINT', () => {
  console.log('\n🛑 Auto-deployment watcher stopped');
  process.exit(0);
}); 