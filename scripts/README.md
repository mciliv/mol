# Deployment Scripts

## Overview
These scripts replace the long command strings in `package.json` with clear, maintainable code.

## Scripts

### `deploy.js` - Production Server Deployment
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

### `ship.js` - Complete Workflow
**Command:** `npm run ship`

**What it does:**
1. ✅ Stages all changes (`git add .`)
2. ✅ Commits with timestamp (`git commit`)
3. ✅ Runs pre-deployment tests
4. ✅ Deploys to Google Cloud Functions
5. ✅ Pushes to git repository (`git push`)

**Same configuration as `deploy.js`**

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