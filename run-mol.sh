#!/bin/bash

# Simple Mol Runner
# ================

echo "🚀 Starting Mol..."

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

# Run tests
echo "🧪 Running tests..."
npm test

# Check if port is already in use
if lsof -i :8080 > /dev/null 2>&1; then
    echo "⚠️  Port 8080 is already in use"
    echo "💡 Killing existing process..."
    pkill -f "node.*server.js" 2>/dev/null || true
    sleep 2
fi

# Start the server
echo "🌐 Starting server..."
echo "📱 The app will automatically find an available port if 8080 is busy"
echo ""
echo "Press Ctrl+C to stop"
echo ""

npm start 