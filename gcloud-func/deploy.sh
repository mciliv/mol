#!/bin/bash

# GCloud Function Deployment Script
# =================================

# Load utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/utils.sh"
source "$SCRIPT_DIR/config.sh"

# Setup logging
setup_logging "$(basename "$0")"

log_info "Starting GCloud Function deployment for $PROJECT_NAME"

# Validate environment
validate_env "PROJECT_NAME" "FUNCTION_NAME" "REGION" "SOURCE_DIR" || {
    log_error "Missing required configuration. Please update config.sh"
    exit 1
}

# Check authentication
check_auth || {
    log_error "Please authenticate with Google Cloud: gcloud auth login"
    exit 1
}

# Deploy function
if [ -d "$SOURCE_DIR" ]; then
    deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR" "$PROJECT_ID" "$RUNTIME" "$MEMORY" "$TIMEOUT" "$ENTRY_POINT"
else
    log_error "Source directory not found: $SOURCE_DIR"
    exit 1
fi

# Test deployment
if [ $? -eq 0 ]; then
    local url=$(get_function_url "$FUNCTION_NAME" "$REGION")
    if [ -n "$url" ]; then
        log_info "Testing deployed function..."
        test_url "$url"
    fi
fi

log_info "GCloud Function deployment completed!" 