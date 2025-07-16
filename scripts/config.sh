#!/bin/bash

# Project Configuration
# ====================

# Domain Configuration
DOMAIN_NAME="queb.space"
DNS_ZONE_NAME="queb-space-zone"
REGION="us-central1"
FUNCTION_NAME="molecular-analysis"

# Google Cloud Configuration
PROJECT_ID="${GOOGLE_CLOUD_PROJECT:-$(gcloud config get-value project 2>/dev/null)}"

# Python Configuration
PYTHON_VERSION="3.12.9"

# Paths
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Export for use in other scripts
export DOMAIN_NAME DNS_ZONE_NAME REGION FUNCTION_NAME PROJECT_ID PYTHON_VERSION SCRIPT_DIR PROJECT_ROOT 