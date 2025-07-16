# Molecular Reality (mol)

## Project Context (for AI Assistants)
**Type**: Molecular visualization web app  
**Key Features**: Camera/photo molecular analysis using OpenAI Vision API  
**Visualization**: Sphere-only representation (NO ball-and-stick), van der Waals radii at 0.8 scale  
**UI Philosophy**: "Keep stupid simple" - minimal design, no extraneous elements  
**Workflow**: Always git commit after acceptance  

This project includes a Node.js server that hosts a static front-end (HTML, CSS, JavaScript) and provides back-end API endpoints for processing molecular data using the OpenAI API, as well as Python utilities for SDF file generation.

## Prerequisites

- [Node.js](https://nodejs.org/) (>= 18)
- [Python](https://www.python.org/) (>= 3.8)
- [poetry](https://python-poetry.org/) (for Python dependencies)

## Python Environment Setup

To set up the Python virtual environment and install dependencies:

```bash
chmod +x scripts/helper.sh
./scripts/helper.sh
```

## Node.js Setup

Install Node.js dependencies (including development dependencies):

```bash
npm install
```

### Environment Variables

Set your OpenAI API key:

```bash
export OPENAI_API_KEY=<your_api_key_here>
```

### Running the Server

Start the Node.js server (serves both front-end and back-end API):

```bash
npm start
```

For development with automatic reloads:

```bash
npm run dev
```

The application will be available at [http://localhost:3000](http://localhost:3000).

## Mobile Setup (iPhone/iPad)

### 1. Start the Server
```bash
npm run dev
```

### 2. Get Your Computer's IP Address
```bash
npm run mobile
```
This will show URLs like:
- HTTPS: `https://192.168.1.100:3001`
- HTTP: `http://192.168.1.100:3000`

### 3. Access on iPhone
1. **Use the HTTPS URL** (iPhone requires HTTPS for camera access)
2. Navigate to `https://[your-ip]:3001` in Safari
3. **Accept the certificate warning**:
   - Tap "Advanced" 
   - Tap "Proceed to [your-ip] (unsafe)"
   - OR tap the "AA" icon in address bar and select "Allow"

### 4. Grant Camera Permission
- When prompted, tap "Allow" for camera access
- If denied, go to Settings > Safari > Camera > Allow

### 5. Troubleshooting
- **Black screen**: Try refreshing the page
- **Permission denied**: Check Safari settings or use the "Try Again" button
- **Certificate issues**: Make sure you're using HTTPS and accepted the certificate
- **Still not working**: Try using ngrok: `npm install -g ngrok && ngrok http 3000`

## Development

### Quick Start
```bash
npm install
npm run dev
```

### Available Scripts
- `npm run dev` - Start development server with hot reload
- `npm run mobile` - Show mobile access URLs
- `npm run cert` - Regenerate SSL certificates
- `npm run tunnel` - Instructions for ngrok tunneling
- `npm run ship` - Complete workflow: commit, test, deploy to production, push to git
- `npm run deploy` - Deploy to Google Cloud Functions (production server only)

### Features
- 📷 Camera-based molecule analysis
- ✏️ Text-based molecule lookup
- 🔍 SMILES string generation
- 🔒 HTTPS support for mobile devices
- 🔄 Hot reload during development

### Technology Stack
- **Frontend**: Vanilla JavaScript, HTML5, CSS3
- **Backend**: Node.js, Express
- **AI**: OpenAI GPT-4 Vision
- **Chemistry**: SMILES notation, molecular analysis
- **Deployment**: Google Cloud Functions

## API Endpoints

- `POST /list-molecules`: Analyze an image click and return possible molecule structures.
- `POST /list-molecules-text`: Provide an object description and get molecule SMILES in JSON.

## Environment Variables
```bash
OPENAI_API_KEY=your_openai_api_key_here
NODE_ENV=development
```

## Deployment

### Prerequisites
- Google Cloud CLI installed and authenticated
- OpenAI API key set in environment

### Deploy to Production

**Complete Workflow (commit, test, deploy to production, push to git):**
```bash
npm run ship
```

**Production Deploy Only (no git changes):**
```bash
npm run deploy
```

The `ship` command:
1. Stages all changes (`git add .`)
2. Commits with timestamp (`git commit`)
3. Runs all tests to ensure code quality
4. **Deploys to Google Cloud Functions (production server)**
5. **Pushes changes to git repository**

### Deployment Structure

The deployment uses `.gcloudignore` to include only essential files:

**Core Application:**
- `index.js` - Cloud Function entry point
- `server.js` - Express server with API endpoints
- `app.js` - Frontend JavaScript application
- `index.html` - Main HTML page
- `style.css` - Application styles
- `schemas.js` - Input validation schemas

**Molecular Processing:**
- `molecular-processor.js` - SMILES processing logic
- `AtomPredictor.js` - AI prediction integration
- `sdf.py`, `crystal.py`, `uniprot.py`, `usdz.py` - Core Python processing modules

**Excluded Development Scripts:**
- `scripts.py` - Development script runner
- `kmeans.py` - Data analysis script
- `jobs.html` - Development interface
- `https-server.js` - Development server
- `lambda.js` - AWS Lambda handler
- `nodemon.json` - Development configuration

**Assets:**
- SVG icons for camera, water, waves

**Benefits:**
- **Faster deployments** - Smaller upload size (~200KB vs full project)
- **Cleaner environment** - No development files or documentation
- **Reduced costs** - Less storage and processing overhead
- **Better security** - No sensitive development files in production
- **Simpler structure** - Flat file organization, no copying needed



## License
MIT
