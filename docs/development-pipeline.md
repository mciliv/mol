# Development Pipeline Documentation

## Infrastructure Overview

### LiveReload System
- **Port**: 35730
- **Middleware**: connect-livereload  
- **Watching**: `backend/**/*`, `frontend/**/*`
- **Extensions**: `.js`, `.html`, `.css`

### Server Configuration
- **HTTP**: `http://localhost:8080`
- **HTTPS**: `https://localhost:3001`
- **Mobile HTTP**: `http://172.20.10.4:8080`
- **Mobile HTTPS**: `https://172.20.10.4:3001`

### Dependencies
**Core**: openai, express, cors, multer, pg, stripe, sharp, fs-extra, zod
**Dev**: nodemon, livereload, connect-livereload, jest

### Database (Optional)
- **Type**: PostgreSQL
- **Name**: mol_users
- **Setup**: `createdb mol_users`

## UI Architecture

### Design Philosophy
- **Style**: "stupid simple" minimalism
- **Rules**: ui.mdc compliance
- **Theme**: Dark (#000000 backgrounds)
- **Interface**: Icons only, no text
- **Elements**: No borders, outlines, box-shadows

### Development Workflow
1. **File Watching**: nodemon restarts server on backend changes
2. **Browser Refresh**: livereload refreshes browser on frontend changes
3. **Auto Debug**: VS Code auto-attach on --inspect
4. **Payment Bypass**: developer mode auto-enabled on localhost
5. **Error Reporting**: structured console logs with 🚨 prefix

## Molecular Analysis
- **AI Provider**: OpenAI Vision API
- **Visualization**: 3DMol.js with sphere representation
- **Scale**: van der Waals radii at 0.8
- **Chemistry Formats**: SMILES, SDF
- **Analysis Types**: text_input, camera_capture, image_upload

## Startup Sequence
1. Load dependencies (openai, express, cors, etc.)
2. Initialize PostgreSQL connection (optional)
3. Start LiveReload server on port 35730
4. Launch HTTP server on port 8080
5. Launch HTTPS server on port 3001
6. Initialize frontend logger and debug systems
7. Auto-enable developer mode for localhost
8. Begin file watching for auto-reload

## Available Commands
- `./dev` - Standard development with live reload
- `./dev-debug` - Maximum error visibility
- `npm test` - Run test suite
- `npm run clean` - Clean up processes 