#!/bin/bash

# Simple local run script for Mol
# ===============================

echo "🚀 Starting Mol locally..."

# Check if Node.js is installed
if ! command -v node &> /dev/null; then
    echo "❌ Node.js not found. Please install Node.js first."
    exit 1
fi

# Check if dependencies are installed
if [ ! -d "node_modules" ]; then
    echo "📦 Installing dependencies..."
    npm install
fi

# Run tests first
echo "🧪 Running tests..."
npm test

# Start the server
echo "🌐 Starting server on http://localhost:8080"
echo "📱 Mobile access: http://$(ifconfig | grep 'inet ' | grep -v 127.0.0.1 | awk '{print $2}' | head -1):8080"
echo ""
echo "Press Ctrl+C to stop"
echo ""

npm start 