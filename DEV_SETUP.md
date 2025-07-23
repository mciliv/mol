# Development Setup Guide

## Quick Start

### Option 1: Simple Development (Recommended)
```bash
./dev
```
This starts the development server with LiveReload enabled.

### Option 2: Without LiveReload (if you're having port conflicts)
```bash
./dev simple
```
This starts the server without LiveReload to avoid port conflicts.

### Option 3: Clean Start (if you have stuck processes)
```bash
./dev clean
```
This kills any existing processes and starts fresh.

## Available Commands

### Development Scripts
- `./dev` - Start with LiveReload (default)
- `./dev simple` - Start without LiveReload
- `./dev watch` - Start with nodemon (backend file watching)
- `./dev clean` - Clean processes and start fresh
- `./dev no-lr` - Start without LiveReload
- `./dev https` - Start with HTTPS
- `./dev setup` - Setup development environment
- `./dev help` - Show help

### NPM Scripts
- `npm run dev` - Start development server with LiveReload
- `npm run dev:simple` - Start without LiveReload
- `npm run dev:watch` - Start with nodemon watching backend files
- `npm run dev:clean` - Kill processes and start fresh
- `npm run dev:no-lr` - Start without LiveReload
- `npm run dev:https` - Start with HTTPS
- `npm run clean` - Kill all related processes

## Troubleshooting

### LiveReload Port Conflicts
If you see "LiveReload port 35729 is already in use":

1. **Quick fix**: Use `./dev simple` to start without LiveReload
2. **Clean fix**: Use `./dev clean` to kill all processes and start fresh
3. **Manual fix**: Run `npm run clean` then `./dev`

### Port Conflicts
If port 8080 is already in use:

1. The new development server automatically finds available ports
2. Use `./dev clean` to kill existing processes
3. Check what's using the port: `lsof -i :8080`

### Database Issues
If you see database connection errors:

1. Make sure PostgreSQL is running: `brew services start postgresql`
2. Create the database: `createdb mol_users`
3. The app will work without database in development mode

## Development Features

### Automatic Port Management
- Automatically finds available ports if default ports are in use
- Kills existing processes on startup to avoid conflicts
- Graceful shutdown handling

### LiveReload
- Watches frontend files for changes
- Automatically refreshes browser
- Configurable file extensions and ignore patterns
- Fallback to polling if file watching fails

### Logging
- All logs go to terminal/server instead of browser console
- Configurable log levels via URL parameter: `?log=debug`
- Structured logging for different event types

### HTTPS Support
- Automatic SSL certificate generation
- HTTPS server on port 3001
- Mobile device access via local IP

## Environment Variables

- `PORT` - Main server port (default: 8080)
- `LIVERELOAD_PORT` - LiveReload port (default: 35729)
- `NODE_ENV` - Environment (development/production)
- `DB_HOST` - Database host (default: localhost)
- `DB_NAME` - Database name (default: mol_users)

## File Structure

```
mol/
├── backend/
│   ├── api/
│   │   ├── server.js          # Main server
│   │   └── dev-server.js      # Development server
│   └── services/
├── frontend/
│   ├── core/                  # Main app files
│   ├── components/            # UI components
│   └── assets/               # CSS, images, etc.
├── dev                        # Development script
├── nodemon.json              # Nodemon configuration
└── package.json              # Dependencies and scripts
```

## Tips

1. **Use `./dev` for normal development** - it handles most issues automatically
2. **Use `./dev simple` if LiveReload is causing problems**
3. **Use `./dev clean` if you have stuck processes**
4. **Check the terminal output** - it shows what ports are being used
5. **Use `?log=debug` in the URL** to see detailed logs in the browser 