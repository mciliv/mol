# Mol - Molecular Analysis with GCloud Func

A molecular analysis app with AI-powered chemical identification, deployed using Google Cloud Functions.

## Commands

### Development
```bash
# Simple scripts (recommended)
./dev.sh                    # Run tests + start with nodemon
./test.sh                   # Run all test suites
./deploy.sh                 # Run tests + deploy to GCP
./format.sh                 # Format code with prettier
./ip.sh                     # Get local IP address

# Direct commands
node server.js              # Start server directly
nodemon server.js           # Start with nodemon (no tests)
node --inspect server.js    # Start with Node.js inspector
jest --watch                # Run tests in watch mode
jest --testPathPattern=unit.test.js --verbose --silent  # Run unit tests only
jest --testPathPattern=integration.test.js --verbose    # Run integration tests only
jest --testPathPattern=system.test.js --verbose --detectOpenHandles  # Run system tests only
python -m pytest tests/ -v  # Run Python tests
python -m debugpy --listen 5678 --wait-for-client -m pytest tests/ -v -s  # Run Python tests with debugger

# NPM commands (only where npm is needed)
npm start                   # Start server (npm wrapper)
npm run dev:unsafe          # Start with nodemon (no tests)
npm run test:debug          # Run tests with debugger (requires npm)
npm run test:watch          # Run tests in watch mode (npm wrapper)
npm run test:fixtures       # Run fixture tests (npm wrapper)
npm run pytest:debug        # Run Python tests with debugger (npm wrapper)
```

### Deployment
```bash
# Direct commands
node ship                   # Deploy using ship script
node scripts/auto-deploy.js # Auto-deploy on changes
gcloud functions deploy molecular-analysis --runtime nodejs20 --trigger-http --allow-unauthenticated --memory 1GB --timeout 540s --entry-point molecularAnalysis --source .  # Deploy to GCP
gcloud functions deploy molecular-analysis --gen2 --runtime nodejs20 --trigger-http --allow-unauthenticated --memory 1GB --timeout 540s --region us-central1 --entry-point molecularAnalysis --source . --set-env-vars OPENAI_API_KEY=$OPENAI_API_KEY --quiet  # Deploy with env vars

# NPM commands (only where npm is needed)
npm run deploy:watch        # Auto-deploy on changes (npm wrapper)
npm run deploy:netlify      # Deploy to Netlify (npm wrapper)
```

### Utilities
```bash
# Direct commands
prettier --write .          # Format code with prettier
ifconfig | grep 'inet ' | grep -v 127.0.0.1 | awk '{print $2}' | head -1  # Get local IP address
echo 'Access on mobile:' && echo 'HTTPS: https://'$(ifconfig | grep 'inet ' | grep -v 127.0.0.1 | awk '{print $2}' | head -1)':3001' && echo 'HTTP: http://'$(ifconfig | grep 'inet ' | grep -v 127.0.0.1 | awk '{print $2}' | head -1)':8080'  # Show mobile access URLs
echo 'Generating fresh SSL certificates...' && rm -rf certs && node -e "require('./server.js')"  # Generate SSL certificates
echo 'Install ngrok: npm install -g ngrok' && echo 'Then run: ngrok http 8080'  # Show ngrok tunnel instructions
jest --testPathPattern=unit.test.js --verbose --silent && echo 'Build complete - static files ready'  # Run tests + build static files

# NPM commands (only where npm is needed)
npm run build               # Run tests + build static files (npm wrapper)
```

### Git & Domain
```bash
# Direct commands
./scripts/check-commit.sh   # Check git commit status
git add . && git commit -m 'Auto-commit: $(date)'  # Auto-commit all changes
./scripts/check-domain-status.sh  # Check domain status
echo 'Run after domain verification:' && echo 'gcloud beta run domain-mappings create --service=molecular-analysis --domain=queb.space --region=us-central1' && echo 'gcloud beta run domain-mappings create --service=molecular-analysis --domain=www.queb.space --region=us-central1'  # Show domain setup commands
```

## What the 4 Webkit Versions Were Doing

Previously, there were 4 versions of webkit that were essentially duplicates:

1. **webkit/** - Original complex toolkit with nested structure
2. **webkit-clean/** - "Clean" version (identical to original)
3. **webkit-proper/** - "Proper" version (identical to original)
4. **webkit-temp/** - Minimal version with just `webkit.sh`

These were consolidated into the new **flat `gcloud-func/` structure** for simplicity.

## New GCloud Func Structure

```
gcloud-func/
├── deploy.sh          # Main deployment script
├── utils.sh           # Core utility functions
├── config.sh          # Configuration (mol-specific)
├── auth.sh            # Authentication utilities
├── domain.sh          # Domain management
├── ssl.sh             # SSL certificate utilities
├── monitor.sh         # Monitoring and logging
└── README.md          # Documentation
```

## Usage

### Basic Deployment
```bash
# Deploy function only
./gcloud-func/deploy.sh
```

### Full Deployment with Domain
```bash
# Deploy with domain, SSL, and monitoring
./deploy-mol.sh
```

### Check Status
```bash
# Check function status
source gcloud-func/utils.sh
source gcloud-func/config.sh
check_function "$FUNCTION_NAME" "$REGION"
```

### Monitor Function
```bash
# Check logs and metrics
source gcloud-func/monitor.sh
check_logs "$FUNCTION_NAME"
get_metrics "$FUNCTION_NAME"
```

## Configuration

The `gcloud-func/config.sh` is pre-configured for the mol project:

- **Function**: `molecular-analysis`
- **Region**: `us-central1`
- **Domain**: `queb.space`
- **Runtime**: `nodejs20`

## Benefits of New Structure

- **Flat**: Easy to copy and use
- **Specific**: Focused on Google Cloud Functions
- **Modular**: Load only what you need
- **No Duplicates**: Single source of truth
- **Well-documented**: Clear examples

## Scripts

- `run-mol.sh` - Run mol (check status and open in browser)
- `deploy-mol.sh` - Full deployment with domain and SSL
- `gcloud-func/deploy.sh` - Basic function deployment
- `gcloud-func/monitor.sh` - Monitoring utilities

## Next Steps

1. **Authenticate**: `gcloud auth login`
2. **Deploy**: `./deploy-mol.sh`
3. **Run**: `./run-mol.sh`
