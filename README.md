# Mol - Molecular Analysis with GCloud Func

A molecular analysis app with AI-powered chemical identification, deployed using Google Cloud Functions.

## Quick Setup

```bash
./setup.sh
source ~/.zshrc  # or ~/.bashrc
```

## Commands

After setup, use these aliases from anywhere:

```bash
dev              # Run tests + start with nodemon
test             # Run all test suites
deploy           # Run tests + deploy to GCP
format           # Format code
```

### All Commands

```bash
# Development
dev, start       # Run tests + start with nodemon
server           # Start server directly
debug            # Start with Node.js inspector
unsafe           # Start with nodemon (no tests)

# Testing
test             # Run all test suites
unit             # Run unit tests only
integration      # Run integration tests only
system           # Run system tests only
watch            # Run tests in watch mode
pytest           # Run Python tests
pytest:debug     # Run Python tests with debugger
fixtures         # Run fixture tests

# Deployment
deploy           # Run tests + deploy to GCP
deploy:now       # Deploy with env vars
ship             # Deploy using ship script
auto-deploy      # Auto-deploy on changes

# Utilities
format           # Format code with prettier
ip               # Get local IP address
mobile           # Show mobile access URLs
cert             # Generate SSL certificates
tunnel           # Show ngrok tunnel instructions
build            # Run tests + build static files

# Git & Domain
check-commit     # Check git commit status
commit           # Auto-commit all changes
domain:status    # Check domain status
domain:setup     # Show domain setup commands
```

## What the 4 Webkit Versions Were Doing

Previously, there were 4 versions of webkit that were essentially duplicates:

1. **webkit/** - Original complex toolkit with nested structure
2. **webkit-clean/** - "Clean" version (identical to original)
3. **webkit-proper/** - "Proper" version (identical to original)
4. **webkit-temp/** - Minimal version with just `webkit.sh`

These were consolidated into the new **flat `gcloud-func/` structure** for simplicity.

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
