#!/usr/bin/env node

// Development server with improved LiveReload
const express = require('express');
const cors = require('cors');
const path = require('path');
const fs = require('fs');
const { spawn } = require('child_process');

// Configuration
const config = {
  port: process.env.PORT || 8080,
  livereloadPort: process.env.LIVERELOAD_PORT || 35729,
  frontendPath: path.join(__dirname, '..', '..', 'frontend'),
  backendPath: path.join(__dirname, '..'),
  enableLiveReload: process.argv.includes('--livereload') !== false,
  enableHttps: process.argv.includes('--https') !== false
};

// Utility functions
const findAvailablePort = async (startPort) => {
  const net = require('net');
  
  const isPortAvailable = (port) => {
    return new Promise((resolve) => {
      const server = net.createServer();
      server.listen(port, () => {
        server.once('close', () => resolve(true));
        server.close();
      });
      server.on('error', () => resolve(false));
    });
  };

  for (let port = startPort; port < startPort + 100; port++) {
    if (await isPortAvailable(port)) {
      return port;
    }
  }
  throw new Error(`No available ports found starting from ${startPort}`);
};

const killProcessOnPort = async (port) => {
  try {
    const { exec } = require('child_process');
    const util = require('util');
    const execAsync = util.promisify(exec);
    
    // Find process using the port
    const { stdout } = await execAsync(`lsof -ti:${port}`);
    if (stdout.trim()) {
      const pids = stdout.trim().split('\n');
      for (const pid of pids) {
        console.log(`🔄 Killing process ${pid} on port ${port}`);
        await execAsync(`kill -9 ${pid}`);
      }
      // Wait a moment for the port to be released
      await new Promise(resolve => setTimeout(resolve, 1000));
    }
  } catch (error) {
    // Process not found or already killed
  }
};

// Start LiveReload server
const startLiveReload = async () => {
  if (!config.enableLiveReload) {
    console.log('⚠️ LiveReload disabled');
    return null;
  }

  try {
    // Try to kill any existing LiveReload process
    await killProcessOnPort(config.livereloadPort);
    
    // Find available port
    const livereloadPort = await findAvailablePort(config.livereloadPort);
    
    console.log(`🔄 Starting LiveReload on port ${livereloadPort}`);
    
    // Start LiveReload server
    const livereload = require('livereload');
    const connectLivereload = require('connect-livereload');
    
    const lrServer = livereload.createServer({
      exts: ['html', 'css', 'js', 'json'],
      ignore: [
        'node_modules/**',
        '**/*.test.js',
        '**/test/**',
        '**/.git/**',
        '**/uploads/**',
        '**/temp/**',
        '**/*.log'
      ],
      port: livereloadPort,
      delay: 100,
      usePolling: true
    });

    // Watch frontend directory
    lrServer.watch(config.frontendPath);
    
    // Error handling
    lrServer.server.on('error', (err) => {
      console.error('❌ LiveReload error:', err.message);
    });

    lrServer.server.on('listening', () => {
      console.log(`✅ LiveReload server running on port ${livereloadPort}`);
      console.log(`👀 Watching: ${config.frontendPath}`);
    });

    return { server: lrServer, port: livereloadPort };
  } catch (error) {
    console.error('❌ Failed to start LiveReload:', error.message);
    return null;
  }
};

// Start main development server
const startDevServer = async () => {
  try {
    // Kill any existing process on the main port
    await killProcessOnPort(config.port);
    
    // Find available port
    const port = await findAvailablePort(config.port);
    
    console.log(`🚀 Starting development server on port ${port}`);
    
    // Start the main server
    const serverProcess = spawn('node', ['backend/api/server.js'], {
      stdio: 'inherit',
      env: {
        ...process.env,
        PORT: port,
        NODE_ENV: 'development',
        LIVERELOAD_PORT: config.livereloadPort
      }
    });

    // Handle server process events
    serverProcess.on('error', (error) => {
      console.error('❌ Server process error:', error);
    });

    serverProcess.on('exit', (code) => {
      console.log(`📴 Server process exited with code ${code}`);
      process.exit(code);
    });

    // Handle graceful shutdown
    const gracefulShutdown = (signal) => {
      console.log(`\n🛑 Received ${signal}, shutting down gracefully...`);
      serverProcess.kill('SIGTERM');
      
      setTimeout(() => {
        console.log('⚠️ Force killing server process...');
        serverProcess.kill('SIGKILL');
        process.exit(1);
      }, 5000);
    };

    process.on('SIGINT', () => gracefulShutdown('SIGINT'));
    process.on('SIGTERM', () => gracefulShutdown('SIGTERM'));

    return { process: serverProcess, port };
  } catch (error) {
    console.error('❌ Failed to start development server:', error.message);
    throw error;
  }
};

// Main function
const main = async () => {
  console.log('🔧 Starting development environment...');
  console.log(`📁 Frontend path: ${config.frontendPath}`);
  console.log(`📁 Backend path: ${config.backendPath}`);
  console.log(`🔄 LiveReload: ${config.enableLiveReload ? 'enabled' : 'disabled'}`);
  console.log(`🔒 HTTPS: ${config.enableHttps ? 'enabled' : 'disabled'}`);

  try {
    // Start LiveReload first
    const livereload = await startLiveReload();
    
    // Start main server
    const server = await startDevServer();
    
    console.log('\n✅ Development environment started successfully!');
    console.log(`🌐 Main server: http://localhost:${server.port}`);
    if (livereload) {
      console.log(`🔄 LiveReload: http://localhost:${livereload.port}`);
    }
    console.log('\n💡 Press Ctrl+C to stop the development server');
    
  } catch (error) {
    console.error('❌ Failed to start development environment:', error.message);
    process.exit(1);
  }
};

// Run if this file is executed directly
if (require.main === module) {
  main();
}

module.exports = { main, config }; 