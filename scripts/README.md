# Scripts Directory

## Overview
This directory contains utility scripts for development, deployment, and infrastructure management. Scripts are designed to be reusable and configurable through a centralized configuration system.

## Configuration Management

### `config.sh` - Centralized Configuration
All scripts use a shared configuration file that centralizes project-specific settings:

- **Domain settings**: Domain name, DNS zone, region
- **Google Cloud settings**: Project ID, function names
- **Python settings**: Version requirements
- **Paths**: Script and project root directories

### Benefits of Centralized Config:
- **🔧 Reusable**: Scripts work across different projects
- **📝 Maintainable**: Single place to update settings
- **🚀 Portable**: Easy to adapt for new projects
- **🐛 Debuggable**: Clear separation of logic and configuration

## Best Practices for Context Management

### 1. **Configuration-First Approach**
- All project-specific values go in `config.sh`
- Scripts should never hardcode domain names, regions, etc.
- Use environment variables with sensible defaults

### 2. **Generic Functions**
- Create reusable functions for common operations
- See `template-generic.sh` for examples
- Functions should accept parameters rather than use global variables

### 3. **Logging and Error Handling**
- All scripts should log their actions
- Use consistent error handling patterns
- Include timestamps and context in logs

### 4. **Documentation**
- Each script should have a clear purpose statement
- Document required environment variables
- Include usage examples

### 5. **Testing**
- Scripts should be testable in isolation
- Use dry-run modes where possible
- Validate configuration before execution

## Scripts

### Infrastructure Scripts

#### `check-domain-status.sh` - Domain Health Check
**Purpose:** Comprehensive domain and DNS status check
**Usage:** `./scripts/check-domain-status.sh`

**What it checks:**
- DNS zone configuration
- Domain verification status
- Cloud Function deployment
- DNS propagation
- Nameserver configuration

**Configuration:** Uses `config.sh` for domain and function settings

#### `helper.sh` - Python Environment Setup
**Purpose:** Sets up Python environment with pyenv and poetry
**Usage:** `source scripts/helper.sh`

**What it does:**
- Installs pyenv if needed
- Sets up Python version from config
- Installs poetry
- Configures virtual environment

### Deployment Scripts

#### `deploy` - Production Server Deployment
**Command:** `npm run deploy`

**What it does:**
1. ✅ Validates environment variables
2. ✅ Runs pre-deployment tests
3. ✅ Deploys to Google Cloud Functions
4. ✅ Shows deployment URL

**Configuration:**
- Function: `molecular-analysis`
- Runtime: `nodejs20`
- Region: `us-central1`
- Memory: `1GB`
- Timeout: `540s`

#### `ship` - Complete Workflow
**Command:** `npm run ship`

**What it does:**
1. ✅ Stages all changes (`git add .`)
2. ✅ Commits with timestamp (`git commit`)
3. ✅ Runs pre-deployment tests
4. ✅ Deploys to Google Cloud Functions
5. ✅ Pushes to git repository (`git push`)

**Same configuration as `deploy`**

## Benefits

- **📖 Readable** - Clear, documented code instead of long strings
- **🔧 Maintainable** - Easy to modify configuration
- **🐛 Debuggable** - Better error handling and logging
- **⚙️ Configurable** - Centralized configuration object
- **📝 Documented** - Clear comments and step-by-step process

## Environment Variables

**Required:**
- `OPENAI_API_KEY` - OpenAI API key for molecular analysis

**Optional:**
- `GOOGLE_CLOUD_PROJECT` - Google Cloud project ID (for URL display)

## Adapting for New Projects

### Quick Setup
1. **Copy the scripts directory** to your new project
2. **Update `config.sh`** with your project-specific values:
   ```bash
   DOMAIN_NAME="your-domain.com"
   DNS_ZONE_NAME="your-dns-zone"
   FUNCTION_NAME="your-function-name"
   ```
3. **Update environment variables** in your deployment
4. **Test scripts** to ensure they work with your configuration

### Customization Points
- **Domain scripts**: Update DNS zone names and domain mappings
- **Function scripts**: Change function names and regions
- **Python scripts**: Adjust version requirements
- **Logging**: Modify log file paths and formats

### Validation
Run `./scripts/check-domain-status.sh` to verify your configuration is working correctly. 