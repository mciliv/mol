# Debug Launcher - VS Code Equivalent

This project provides functional equivalents to VS Code's `launch.json` and `tasks.json` configurations as standalone scripts. You can use these scripts to debug your full-stack application without needing VS Code.

## Files

- `debug-launcher.sh` - Bash version of the debug launcher
- `debug_launcher.py` - Python version of the debug launcher  
- `debug` - Wrapper script that can run either version
- `DEBUG_LAUNCHER.md` - This documentation

## Quick Start

```bash
# Make scripts executable (if not already)
chmod +x debug debug-launcher.sh debug_launcher.py

# Run interactive menu (bash version)
./debug

# Run interactive menu (Python version)
./debug --python

# Run specific mode directly
./debug quick-debug
./debug --python full-stack-hot-reload
```

## Available Modes

### Compound Configurations (VS Code Compounds)

1. **🚀 Full Stack Hot Reload Debug** - Starts nodemon + Chrome + Python debugging
2. **🔗 Full Stack Attach Debug** - Attaches to existing processes
3. **⚡ Quick Debug (Backend + Frontend)** - Starts nodemon + Chrome debugging
4. **🐍 Python + Node Debug** - Starts nodemon + Python debugging

### Individual Components

5. **🚀 Start Node.js Server** - Direct Node.js server with debugging
6. **🔗 Attach to Nodemon** - Connect to existing nodemon process
7. **🌐 Start Chrome with debugging** - Chrome with remote debugging
8. **🌐 Start Edge with debugging** - Edge with remote debugging
9. **🦊 Start Firefox with debugging** - Firefox with remote debugging
10. **🔗 Attach to Chrome** - Connect to existing Chrome process
11. **🐍 Start Python debugging** - Python debugpy server
12. **🐍 Run Python tests** - Python tests with debugging
13. **🐍 Run Python script** - Custom Python script with debugging
14. **🔗 Attach to Python** - Connect to existing Python process
15. **🧪 Run Jest tests** - Jest tests in watch mode

### Utilities

16. **📊 Show status** - Display running processes and ports
17. **🛑 Stop all processes** - Clean up all debug processes
18. **❌ Exit** - Exit the launcher

## Usage Examples

### Interactive Menu
```bash
./debug                    # Bash version
./debug --python          # Python version
```

### Direct Mode Execution
```bash
# Start full stack development
./debug full-stack-hot-reload

# Quick backend + frontend debug
./debug quick-debug

# Python debugging only
./debug python-debug

# Run Python tests with debugging
./debug python-tests

# Run specific Python script
./debug python-script sdf.py --help

# Check status
./debug status

# Stop all processes
./debug stop
```

### Python Version Specific
```bash
# Interactive menu
./debug --python

# Direct execution
./debug --python full-stack-hot-reload
./debug --python python-script --script uniprot.py --args "--version"
```

## Debug Ports

The launcher uses the following default ports:

- **Node.js Debug**: `localhost:9229`
- **Chrome Debug**: `localhost:9223` 
- **Python Debug**: `localhost:5678`
- **App Server**: `localhost:3000`

## Connecting Debuggers

### Node.js Debugging
Connect your Node.js debugger to `localhost:9229`

### Chrome Debugging
1. Open Chrome DevTools
2. Go to `chrome://inspect`
3. Click "Configure" and add `localhost:9223`
4. Click "inspect" on your app

### Python Debugging
Connect your Python debugger to `localhost:5678`

## Environment Variables

The launcher sets these environment variables:

- `PY_DEBUG=1` - Enables Python debugging
- `NODE_ENV=development` - Sets Node.js environment
- `PYTHONPATH` - Set to workspace root for Python imports

## VS Code Launch.json Equivalents

| VS Code Configuration | Script Equivalent |
|----------------------|-------------------|
| `Node.js Server` | `./debug node-server` |
| `Attach to Nodemon` | `./debug attach-nodemon` |
| `Chrome Debug` | `./debug chrome` |
| `Chrome Attach` | `./debug attach-chrome` |
| `Python Debug` | `./debug python-tests` |
| `Python Attach` | `./debug attach-python` |
| `Python Script Debug` | `./debug python-script sdf.py --help` |
| `Jest Tests` | `./debug jest` |
| `Edge Debug` | `./debug edge` |
| `Firefox Debug` | `./debug firefox` |

## VS Code Compounds Equivalents

| VS Code Compound | Script Equivalent |
|------------------|-------------------|
| `🚀 Full Stack Hot Reload Debug` | `./debug full-stack-hot-reload` |
| `🔗 Full Stack Attach Debug` | `./debug full-stack-attach` |
| `⚡ Quick Debug (Backend + Frontend)` | `./debug quick-debug` |
| `🐍 Python + Node Debug` | `./debug python-node` |

## VS Code Tasks.json Equivalents

| VS Code Task | Script Equivalent |
|--------------|-------------------|
| `start-nodemon` | `npm run dev` (handled by launcher) |
| `delay` | Built into launcher timing |
| `wait-for-server` | Built into launcher port checking |

## Requirements

### For Bash Version
- Bash shell
- Node.js and npm
- nodemon (auto-installed if missing)
- Chrome/Chromium, Edge, or Firefox
- Python with debugpy (auto-installed if missing)

### For Python Version
- Python 3.6+
- Node.js and npm
- nodemon (auto-installed if missing)
- Chrome/Chromium, Edge, or Firefox
- debugpy (auto-installed if missing)

## Troubleshooting

### Port Already in Use
If you get port conflicts, the launcher will show an error. You can:
1. Stop existing processes: `./debug stop`
2. Kill processes manually: `lsof -ti:9229 | xargs kill -9`
3. Change ports in the script variables

### Browser Not Found
The launcher will try to find browsers automatically. If it fails:
1. Install the browser
2. Add the browser path to the script
3. Use a different browser

### Python Debugpy Issues
If Python debugging doesn't work:
1. Install debugpy: `pip install debugpy`
2. Check Python path: `python -c "import debugpy"`
3. Verify port availability

### Node.js Debugging Issues
If Node.js debugging doesn't work:
1. Check Node.js version: `node --version`
2. Verify nodemon installation: `nodemon --version`
3. Check if port 9229 is available

## Advanced Usage

### Custom Ports
Edit the script variables to change default ports:
```bash
# In debug-launcher.sh
NODE_DEBUG_PORT=9229
CHROME_DEBUG_PORT=9223
PYTHON_DEBUG_PORT=5678
```

### Custom Scripts
Run any Python script with debugging:
```bash
./debug python-script my_script.py --arg1 --arg2
```

### Background Operation
The launcher runs processes in the background. Use `./debug status` to check running processes and `./debug stop` to clean up.

## Integration with IDEs

### VS Code
You can still use VS Code's debugger by connecting to the ports opened by the launcher:
- Node.js: `localhost:9229`
- Chrome: `localhost:9223`
- Python: `localhost:5678`

### PyCharm
Configure Python debugger to connect to `localhost:5678`

### WebStorm/IntelliJ
Configure Node.js debugger to connect to `localhost:9229`

## Contributing

To add new debug configurations:
1. Add the function to the launcher script
2. Add the menu option
3. Add the command-line argument
4. Update this documentation

## License

This project is part of the mol workspace and follows the same license terms. 