#!/bin/bash

# Run Mol with GCloud Func
# ========================

# Load gcloud-func utilities
source gcloud-func/utils.sh
source gcloud-func/config.sh

# Setup logging
setup_logging "$(basename "$0")"

log_info "Starting Mol with GCloud Func..."

# Check authentication
check_auth || {
    log_error "Please authenticate with Google Cloud first"
    log_info "Run: gcloud auth login"
    exit 1
}

# Check function status
log_info "Checking function status..."
check_function "$FUNCTION_NAME" "$REGION"

# Get function URL
local url=$(get_function_url "$FUNCTION_NAME" "$REGION")
if [ -n "$url" ]; then
    log_info "Function URL: $url"
    
    # Test the function
    log_info "Testing function..."
    test_url "$url"
    
    # Open in browser (macOS)
    if command -v open &> /dev/null; then
        log_info "Opening in browser..."
        open "$url"
    fi
else
    log_error "Function not deployed. Run: ./gcloud-func/deploy.sh"
    exit 1
fi

log_info "Mol is running with GCloud Func!" 