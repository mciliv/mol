# Mol - Molecular Analysis with GCloud Func

A molecular analysis app with AI-powered chemical identification, deployed using Google Cloud Functions.

## Quick Start

### Run Mol
```bash
./run-mol.sh
```

### Deploy Mol
```bash
./deploy-mol.sh
```

## Commands

### Development
```bash
npm start          # Start server
npm run dev        # Run tests + start with nodemon
npm run dev:unsafe # Start with nodemon (no tests)
npm run debug      # Start with Node.js inspector
```

### Testing
```bash
npm test           # Run unit + integration tests
npm run test:watch # Run tests in watch mode
npm run test:debug # Run tests with debugger
npm run test:unit  # Run unit tests only
npm run test:integration # Run integration tests only
npm run test:system # Run system tests only
npm run test:all   # Run all test suites
npm run test:pre-deploy # Run all tests + pytest
npm run pytest     # Run Python tests
npm run pytest:debug # Run Python tests with debugger
```

### Deployment
```bash
npm run deploy:gcp # Deploy to Google Cloud Functions
npm run deploy:now # Deploy with env vars
npm run deploy     # Deploy using ship script
npm run deploy:watch # Auto-deploy on changes
npm run pre-deploy # Check commit + run tests
```

### Utilities
```bash
npm run format     # Format code with prettier
npm run ip         # Get local IP address
npm run mobile     # Show mobile access URLs
npm run cert       # Generate SSL certificates
npm run tunnel     # Show ngrok tunnel instructions
npm run build      # Run tests + build static files
```

### Git & Domain
```bash
npm run check-commit # Check git commit status
npm run commit-all   # Auto-commit all changes
npm run domain:status # Check domain status
npm run domain:setup  # Show domain setup commands
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
