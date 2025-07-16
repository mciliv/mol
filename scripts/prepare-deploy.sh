#!/bin/bash

# Prepare deployment directory with only essential files
echo "🚀 Preparing deployment directory..."

# Clean deploy directory
rm -rf deploy/*
rm -rf deploy/.* 2>/dev/null || true

# Copy essential application files
echo "📁 Copying core application files..."
cp index.js deploy/
cp server.js deploy/
cp app.js deploy/
cp index.html deploy/
cp style.css deploy/
cp schemas.js deploy/

# Copy molecular processing files
echo "🧬 Copying molecular processing files..."
cp molecular-processor.js deploy/
cp AtomPredictor.js deploy/
cp sdf.py deploy/
cp crystal.py deploy/
cp uniprot.py deploy/
cp usdz.py deploy/

# Note: Excluding development scripts (scripts.py, kmeans.py, etc.)

# Copy assets
echo "🎨 Copying assets..."
cp *.svg deploy/

# Copy configuration files
echo "⚙️  Copying configuration files..."
cp package.json deploy/
cp __init__.py deploy/

# Create sdf_files directory (empty)
mkdir -p deploy/sdf_files

echo "✅ Deployment directory prepared!"
echo "📊 Files in deploy/:"
ls -la deploy/

echo ""
echo "🚀 Ready to deploy with:"
echo "cd deploy && gcloud functions deploy molecular-analysis --gen2 --runtime nodejs20 --trigger-http --allow-unauthenticated --memory 1GB --timeout 540s --region us-central1 --entry-point molecularAnalysis --source . --set-env-vars OPENAI_API_KEY=\$OPENAI_API_KEY --quiet" 