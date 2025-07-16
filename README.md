# WebKit - Shared Development Toolkit

A collection of reusable utilities for Google Cloud projects, designed to be shared across multiple applications like a web engine.

## Overview

WebKit provides standardized utilities for:
- **Google Cloud Functions** deployment and management
- **Domain management** and DNS configuration
- **Logging** with consistent formatting and levels
- **Python environment** setup and management
- **Project configuration** templates

## Quick Start

### 1. Clone WebKit
```bash
git clone <webkit-repo-url> webkit
cd webkit
```

### 2. Use in Your Project
```bash
# Source the utilities
source webkit/webkit.sh

# Use the functions
log_info "Hello from WebKit!"
check_gcloud_auth
```

### 3. Create Project Scripts
Copy templates from `templates/` and customize:
```bash
cp templates/deploy.sh.template myproject/scripts/deploy.sh
cp config/config.sh.template myproject/config.sh
```

## Structure

```
webkit/
├── core/                    # Core utility functions
│   ├── logging-utils.sh    # Logging and error handling
│   ├── gcloud-utils.sh     # Google Cloud operations
│   ├── domain-utils.sh     # Domain and DNS management
│   └── python-utils.sh     # Python environment setup
├── templates/              # Script templates
│   ├── deploy.sh.template  # Deployment script template
│   └── domain-check.sh.template # Domain verification template
├── config/                 # Configuration templates
│   └── config.sh.template  # Project configuration template
├── webkit.sh              # Main loader script
└── README.md              # This file
```

## Core Utilities

### Logging Utilities (`core/logging-utils.sh`)
- `log_info()`, `log_warn()`, `log_error()`, `log_debug()`
- `setup_logging()` - Initialize logging for a script
- `set_error_handling()` - Setup error handling and logging
- `show_progress()` - Display progress bars

### Google Cloud Utilities (`core/gcloud-utils.sh`)
- `check_gcloud_auth()` - Verify authentication
- `deploy_function()` - Deploy Cloud Functions
- `check_cloud_function()` - Check function status
- `get_function_url()` - Get deployed function URL
- `validate_gcloud_env()` - Validate environment variables

### Domain Utilities (`core/domain-utils.sh`)
- `check_domain_resolution()` - Test DNS resolution
- `check_gcloud_dns_zone()` - Verify DNS zone exists
- `check_domain_verification()` - Check domain verification status
- `get_gcloud_nameservers()` - Get nameserver configuration
- `check_domain_mappings()` - Check Cloud Run domain mappings
- `test_url()` - Test URL accessibility

### Python Utilities (`core/python-utils.sh`)
- `check_pyenv()` - Verify pyenv installation
- `install_python_version()` - Install specific Python version
- `check_poetry()` - Verify poetry installation
- `setup_python_env()` - Complete Python environment setup
- `check_python_deps()` - Validate dependencies
- `run_python_tests()` - Run Python tests

## Usage Examples

### Basic Script with Logging
```bash
#!/bin/bash
source webkit/webkit.sh

# Setup error handling
set_error_handling "my-script.sh"

# Load project config
source config.sh

log_info "Starting deployment for $PROJECT_NAME"
# ... your script logic
```

### Domain Check Script
```bash
#!/bin/bash
source webkit/webkit.sh
source config.sh

echo "🌐 Checking domain: $DOMAIN_NAME"
check_domain_resolution "$DOMAIN_NAME"
check_gcloud_dns_zone "$DNS_ZONE_NAME"
check_domain_verification "$DOMAIN_NAME"
```

### Deployment Script
```bash
#!/bin/bash
source webkit/webkit.sh
source config.sh

if check_gcloud_auth; then
    deploy_function "$FUNCTION_NAME" "$REGION" "$SOURCE_DIR"
else
    log_error "Not authenticated"
    exit 1
fi
```

## Configuration

Create a `config.sh` file in your project:

```bash
# Copy the template
cp webkit/config/config.sh.template config.sh

# Edit with your project details
PROJECT_NAME="my-app"
PROJECT_ID="my-gcp-project"
REGION="us-central1"
FUNCTION_NAME="my-function"
DOMAIN_NAME="myapp.com"
```

## Templates

### Deployment Template
The `templates/deploy.sh.template` provides a complete deployment script with:
- Authentication checks
- Environment validation
- Test execution
- Function deployment
- URL testing
- Domain mapping instructions

### Domain Check Template
The `templates/domain-check.sh.template` provides comprehensive domain verification:
- DNS zone status
- Domain verification
- Nameserver configuration
- Propagation testing
- Actionable next steps

## Best Practices

1. **Always source webkit.sh** at the start of your scripts
2. **Use set_error_handling()** for proper error management
3. **Load project config** from a separate config.sh file
4. **Use log functions** instead of echo for consistent output
5. **Validate environment** before critical operations
6. **Test deployments** after completion

## Integration with Projects

### Simple Integration
```bash
# Clone webkit alongside your project
git clone <webkit-repo> webkit

# Source in your scripts
source webkit/webkit.sh
```

### Git Submodule Integration
```bash
# Add as submodule
git submodule add <webkit-repo> webkit

# Update submodule
git submodule update --remote webkit
```

## Contributing

1. Add new utilities to `core/`
2. Create templates in `templates/`
3. Update documentation
4. Test with multiple projects
5. Keep utilities generic and reusable

## Version History

- **1.0.0** - Initial release with core utilities
  - Logging system
  - Google Cloud utilities
  - Domain management
  - Python environment setup
  - Deployment and domain check templates

## License

MIT License - see LICENSE file for details.
