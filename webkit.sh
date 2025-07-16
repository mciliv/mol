#!/bin/bash

# Minimal webkit - just the essentials
# ===================================

# Colors
GREEN="\033[32m"
RED="\033[31m"
RESET="\033[0m"

# Simple logging
log() {
    echo -e "${GREEN}[$(date '+%H:%M:%S')] $1${RESET}"
}

error() {
    echo -e "${RED}[$(date '+%H:%M:%S')] ERROR: $1${RESET}"
}

# Basic gcloud helpers
gcloud_auth() {
    gcloud auth list --filter=status:ACTIVE --format="value(account)" 2>/dev/null | grep -q .
}

gcloud_project() {
    gcloud config get-value project 2>/dev/null
}

gcloud_deploy() {
    local name="$1"
    local region="$2"
    local source="$3"
    
    gcloud functions deploy "$name" \
        --region="$region" \
        --source="$source" \
        --runtime=nodejs20 \
        --trigger-http \
        --allow-unauthenticated
}

# Error handling
set -e
trap 'error "Script failed at line $LINENO"; exit 1' ERR 