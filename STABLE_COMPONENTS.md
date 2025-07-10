# Stable Components - DO NOT MODIFY

## ✅ Core Engine Components (STABLE)

### 1. **SMILES Processing Engine**
- **Files**: `sdf.py`, `molecular-processor.js`
- **Status**: ✅ STABLE - Working perfectly
- **Function**: Converts SMILES notation to 3D SDF files
- **Last Modified**: Current version
- **⚠️ WARNING**: Do not modify these files - they handle all molecular generation

### 2. **Schemas & Data Structure**
- **Files**: `schemas.js`
- **Status**: ✅ STABLE - SMILES-only approach working
- **Function**: Defines API schemas for SMILES arrays
- **Last Modified**: Simplified to SMILES-only
- **⚠️ WARNING**: Schema changes break the entire API

### 3. **HTTPS Server Module**
- **Files**: `https-server.js`
- **Status**: ✅ STABLE - SSL certificates and networking working
- **Function**: Handles HTTPS setup with self-signed certificates
- **Last Modified**: Separated from main server
- **⚠️ WARNING**: SSL changes break mobile camera access

### 4. **AI Analyzer Module**
- **Files**: `ai-analyzer.js`
- **Status**: ✅ STABLE - People handling working correctly
- **Function**: OpenAI Vision API integration with human body fallback
- **Last Modified**: Fixed people detection with generic chemistry
- **⚠️ WARNING**: Prompt changes can break video feed analysis

## 🔧 Configuration Components (MODIFY WITH CARE)

### 5. **Main Server**
- **Files**: `server.js`
- **Status**: ⚠️ MODIFY WITH CARE - Core routes stable
- **Function**: Express server with API routes
- **Safe to modify**: Port settings, new routes
- **⚠️ DANGER ZONES**: Existing routes, middleware order

### 6. **Frontend**
- **Files**: `index.html`, `app.js`, `style.css`
- **Status**: ⚠️ MODIFY WITH CARE - Camera integration delicate
- **Function**: Video feed, clicking, 3D visualization
- **Safe to modify**: Styling, UI text
- **⚠️ DANGER ZONES**: Camera permissions, click handlers

## 🧪 Testing Strategy

### Regression Prevention
1. **Before any changes**: Test these core functions:
   - Video feed clicking (materials + people)
   - SMILES generation from text
   - SDF file creation
   - 3D visualization

2. **Component isolation**: 
   - Each component has single responsibility
   - Changes to one component shouldn't affect others
   - Test interfaces between components

3. **Commit stable states**:
   - Commit after each working feature
   - Tag stable releases
   - Document what's working

## 📋 Change Protocol

### Before modifying ANY component:
1. ✅ Test current functionality
2. ✅ Commit current state
3. ✅ Make minimal changes
4. ✅ Test immediately
5. ✅ Commit if working
6. ❌ If broken, revert immediately

### Emergency Revert Commands:
```bash
# Revert last commit
git reset --hard HEAD~1

# Restart server cleanly
pkill -f "node server.js"
npm run dev
```

## 🎯 Working Features (DO NOT BREAK)

### Core Functionality ✅
- [x] Video feed camera access
- [x] Click to analyze materials
- [x] Click to analyze people → human body chemistry
- [x] Text input analysis
- [x] SMILES → SDF generation
- [x] 3D molecular visualization
- [x] HTTPS for mobile access

### API Endpoints ✅
- [x] `POST /image-molecules` - Video analysis
- [x] `POST /object-molecules` - Text analysis  
- [x] `POST /generate-sdfs` - SMILES array to SDF
- [x] `GET /sdf_files/*` - SDF file serving

## 🚨 Red Flags - Stop Development If You See:

1. **"Failed to parse AI response"** - AI analyzer broken
2. **"EADDRINUSE"** errors - Port conflicts
3. **"SMILES generation failed"** - Molecular engine broken
4. **"Camera not accessible"** - HTTPS/permissions broken
5. **Black 3D viewer screens** - Visualization broken

When you see these: **STOP, COMMIT CURRENT STATE, REVERT TO WORKING VERSION** 