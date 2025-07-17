#!/bin/bash

# Deploy Mol with GCloud Func
# ===========================

# Load gcloud-func utilities
source gcloud-func/utils.sh
source gcloud-func/config.sh
source gcloud-func/auth.sh
source gcloud-func/domain.sh
source gcloud-func/ssl.sh
source gcloud-func/monitor.sh

# Setup logging
setup_logging "$(basename "$0")"

log_info "Deploying Mol with GCloud Func..."

# Verify setup
verify_setup || {
    log_error "Setup verification failed"
    exit 1
}

# Deploy function
log_info "Deploying function: $FUNCTION_NAME"
deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"

if [ $? -eq 0 ]; then
    log_info "✅ Function deployed successfully"
    
    # Setup domain if configured
    if [ -n "$DOMAIN_NAME" ] && [ "$DOMAIN_NAME" != "your-domain.com" ]; then
        log_info "Setting up domain: $DOMAIN_NAME"
        setup_domain "$DOMAIN_NAME"
        
        # Setup SSL
        log_info "Setting up SSL certificate"
        setup_ssl "$DOMAIN_NAME"
        
        # Wait for SSL
        log_info "Waiting for SSL certificate to be active..."
        wait_for_ssl "$DOMAIN_NAME"
        
        # Test domain
        log_info "Testing domain..."
        test_url "https://$DOMAIN_NAME"
    fi
    
    # Setup monitoring
    log_info "Setting up monitoring..."
    setup_monitoring "$FUNCTION_NAME"
    setup_alerts "$FUNCTION_NAME"
    
    # Get function URL
    local url=$(get_function_url "$FUNCTION_NAME" "$REGION")
    log_info "✅ Mol deployed successfully!"
    log_info "Function URL: $url"
    
    if [ -n "$DOMAIN_NAME" ] && [ "$DOMAIN_NAME" != "your-domain.com" ]; then
        log_info "Domain URL: https://$DOMAIN_NAME"
    fi
else
    log_error "❌ Function deployment failed"
    exit 1
fi 