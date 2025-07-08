# 🚀 Full-Stack Hot Reload Debugging Guide

## Quick Start

### 1. **🚀 Full Stack Hot Reload Debug** (Recommended)
- **What it does**: Starts nodemon + Chrome + Python debugging simultaneously
- **Hot reload**: All levels of the stack reload automatically on file changes
- **Usage**: Select this from debug dropdown and press F5

### 2. **🔗 Full Stack Attach Debug**
- **What it does**: Attaches to already running processes
- **Use when**: You have nodemon already running
- **Usage**: Start nodemon manually, then select this configuration

### 3. **⚡ Quick Debug (Backend + Frontend)**
- **What it does**: Just Node.js + Chrome debugging
- **Use when**: You only need frontend/backend debugging
- **Usage**: Fastest option for most debugging scenarios

## Debugging Ports

| Port | Purpose | Process |
|------|---------|---------|
| 9229 | Node.js Backend | Nodemon with `--inspect=0.0.0.0:9229` |
| 9223 | Chrome Frontend | Chrome DevTools for `app.js` |
| 5678 | Python Scripts | Python debugpy for `sdf.py`, etc. |
| 3000 | Your App | `http://localhost:3000` |

## Setting Breakpoints

### Backend (Node.js)
- **Files**: `server.js`, `schemas.js`
- **Good spots**: 
  - Line 58: `/image-molecules` route start
  - Line 95: `/object-molecules` route start
  - Line 85: After OpenAI API call

### Frontend (Chrome)
- **Files**: `app.js`, `index.html`
- **Good spots**:
  - Line 50: `handleInteraction` function
  - Line 35: Text input handler
  - Line 120: Camera setup

### Python Scripts
- **Files**: `sdf.py`, `dock.py`, `chem.py`
- **Good spots**: Main functions, API calls

## Hot Reload Features

### ✅ What Reloads Automatically
- **Backend**: Nodemon watches `*.js`, `*.json`, `*.html`, `*.css`
- **Frontend**: Livereload refreshes browser on file changes
- **Debugger**: Stays connected through restarts

### 🔄 Manual Restart Commands
```bash
# In nodemon terminal
rs                    # Restart nodemon
Ctrl+C, npm run dev   # Full restart

# In Cursor
F5                    # Restart debug session
Ctrl+Shift+F5         # Stop and restart
```

## Debugging Workflow

### 1. Start Full Stack Debug
1. Open Cursor debug panel (`Cmd+Shift+D`)
2. Select "🚀 Full Stack Hot Reload Debug"
3. Press F5 or click green play button
4. Wait for all processes to start

### 2. Set Breakpoints
1. Open `server.js` - set breakpoint on line 58
2. Open `app.js` - set breakpoint on line 50
3. Open `sdf.py` - set breakpoint in main function

### 3. Test the Flow
1. Open `http://localhost:3000` in browser
2. Click on video feed or enter text
3. Debugger will pause at your breakpoints
4. Step through code with F10 (step over), F11 (step into)

### 4. Hot Reload Testing
1. Make a change to `server.js`
2. Save the file
3. Nodemon restarts automatically
4. Debugger reconnects automatically
5. Test your changes immediately

## Troubleshooting

### Breakpoints Not Hitting
1. **Check debugger is attached**: Look for "Attach to Nodemon" in debug panel
2. **Verify ports**: `lsof -i :9229` should show Node.js process
3. **Restart debugger**: Stop and restart the debug session

### Hot Reload Not Working
1. **Check nodemon**: Should show "restarting due to changes"
2. **Check livereload**: Browser should refresh automatically
3. **File watching**: Ensure files are in watched directories

### Port Conflicts
1. **Kill existing processes**: `pkill -f nodemon`
2. **Check ports**: `lsof -i :9229` and `lsof -i :9223`
3. **Restart everything**: Stop all debug sessions and restart

## Advanced Features

### Conditional Breakpoints
- Right-click on breakpoint → Edit breakpoint
- Add condition: `x > 100` or `object === 'aspirin'`

### Logpoints
- Right-click on breakpoint → Add logpoint
- Log without stopping: `Click at (${x}, ${y})`

### Watch Expressions
- Add variables to watch in debug panel
- Monitor `req.body`, `parsed.output_parsed`, etc.

### Call Stack Navigation
- Use call stack to navigate through function calls
- See the full execution path when paused

## Performance Tips

### Fast Debugging
1. Use "⚡ Quick Debug" for most scenarios
2. Set breakpoints only where needed
3. Use logpoints instead of breakpoints for logging

### Memory Management
1. Stop debug sessions when not needed
2. Clear browser cache if frontend debugging gets slow
3. Restart nodemon if backend gets sluggish

## File Structure for Debugging

```
mol/
├── server.js          # Backend API routes (port 9229)
├── app.js            # Frontend logic (port 9223)
├── sdf.py            # Python scripts (port 5678)
├── schemas.js        # Zod schemas
├── .vscode/
│   ├── launch.json   # Debug configurations
│   └── tasks.json    # Pre-launch tasks
└── DEBUGGING.md      # This guide
```

## Quick Commands

```bash
# Start nodemon manually
npm run dev

# Check debug ports
lsof -i :9229  # Node.js
lsof -i :9223  # Chrome
lsof -i :5678  # Python

# Kill all debug processes
pkill -f nodemon
pkill -f chrome.*remote-debugging
```

Happy debugging! 🐛✨ 