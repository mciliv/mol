# Mol - Molecular Analysis App

A simplified molecular visualization app with AI-powered chemical identification.

## Quick Start

```bash
# Install dependencies
npm install

# Start development
./run dev

# Run tests
./run test

# Deploy to production
./run deploy
```

## Commands

The `./run` script provides all functionality:

```bash
# Development
./run dev          # Run tests + start development server
./run server       # Start production server
./run debug        # Start server with debugger
./run cleanup      # Clean up processes and ports

# Testing
./run test         # Run all tests (JavaScript + Python)
./run test:unit    # Run JavaScript tests only
./run test:python  # Run Python tests only

# Deployment
./run deploy       # Deploy to Google Cloud Functions
./run ship         # Complete workflow: commit, test, deploy, push
```

## Environment Variables

For deployment:
- `OPENAI_API_KEY` - Required for AI molecular analysis

## Architecture

- **Frontend**: Vanilla JavaScript with 3Dmol.js for visualization
- **Backend**: Node.js/Express API
- **AI**: OpenAI Vision API for molecular analysis
- **Chemistry**: Python scripts for SDF processing
- **Deployment**: Google Cloud Functions

## Development

Local development auto-enables developer mode. For production, users need to set up payment.

The app uses sphere representation only for molecules (van der Waals radii at 0.8 scale).

## Simplified Structure

This codebase has been streamlined by:
- Consolidating multiple run scripts into one
- Removing duplicate package.json files
- Cleaning up deprecated payment methods
- Removing complex infrastructure tooling
- Simplifying the VS Code debug configuration

Core functionality remains unchanged while reducing complexity.