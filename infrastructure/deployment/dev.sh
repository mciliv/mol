#!/bin/bash

# dev.sh - Development server startup with comprehensive setup

set -e

# Get the absolute path to project root
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_ROOT"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}ℹ️  $1${NC}"
}

log_success() {
    echo -e "${GREEN}✅ $1${NC}"
}

log_warning() {
    echo -e "${YELLOW}⚠️  $1${NC}"
}

# Check if .env file exists and database is set up
check_database_setup() {
    if [ ! -f ".env" ]; then
        log_warning "No .env file found. Setting up database..."
        ./infrastructure/scripts/setup-database.sh
    else
        # Source .env file
        export $(grep -v '^#' .env | xargs)
        
        # Check if database connection works
        if command -v psql >/dev/null 2>&1; then
            DB_HOST=${DB_HOST:-localhost}
            DB_PORT=${DB_PORT:-5432}
            DB_NAME=${DB_NAME:-mol_users}
            DB_USER=${DB_USER:-mol_user}
            DB_PASSWORD=${DB_PASSWORD:-mol_password}
            
            if ! PGPASSWORD=$DB_PASSWORD psql -h $DB_HOST -p $DB_PORT -U $DB_USER -d $DB_NAME -c 'SELECT 1;' >/dev/null 2>&1; then
                log_warning "Database connection failed. Running setup..."
                ./infrastructure/scripts/setup-database.sh
            else
                log_success "Database connection verified"
            fi
        fi
    fi
}

echo "🧹 Cleaning up any existing development processes..."
./cleanup 2>/dev/null || true

echo "🔧 Checking database setup..."
check_database_setup

echo "🚀 Starting development server..."

# Load environment variables if .env exists
if [ -f ".env" ]; then
    export $(grep -v '^#' .env | xargs)
fi

# Run tests first to ensure everything is working
npm test --prefix infrastructure/config

# Start the development server with nodemon
cd infrastructure/config
npx nodemon ../../backend/api/server.js 