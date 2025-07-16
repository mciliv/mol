# Minimal Webkit

Just the essentials for Google Cloud projects.

## Usage

```bash
# Source the webkit
source webkit.sh

# Use the functions
log "Starting deployment..."
gcloud_deploy "my-function" "us-central1" "./src"
```

## Functions

- `log "message"` - Simple green logging
- `error "message"` - Red error logging  
- `gcloud_auth` - Check if authenticated
- `gcloud_project` - Get current project
- `gcloud_deploy name region source` - Deploy function

## Git Submodule

```bash
# Add as submodule
git submodule add https://github.com/your-org/webkit.git

# Use in scripts
source webkit/webkit.sh
```

That's it. No bureaucracy. 