#!/bin/bash

# Debug Launcher Script - Functional equivalent to VS Code launch.json & tasks.json
# Usage: ./debug-launcher.sh [mode] [options]

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
PURPLE='\033[0;35m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Configuration
WORKSPACE_ROOT="$(pwd)"
NODE_PORT=${PORT:-8080}  # Use environment PORT or default to 8080 (Cloud Functions standard)
NODE_DEBUG_PORT=9229
CHROME_DEBUG_PORT=9223
PYTHON_DEBUG_PORT=5678

# Process management (compatible with older bash versions)
PIDS=""
PROCESS_NAMES=""

# Cleanup function
cleanup() {
    echo -e "${YELLOW}🛑 Cleaning up processes...${NC}"
    if [ -n "$PIDS" ]; then
        for pid in $PIDS; do
            if kill -0 "$pid" 2>/dev/null; then
                echo -e "${YELLOW}Killing process $pid${NC}"
                kill "$pid" 2>/dev/null || true
            fi
        done
    fi
    exit 0
}

# Trap cleanup on script exit
trap cleanup EXIT INT TERM

# Helper functions
log() {
    echo -e "${GREEN}[$(date +'%H:%M:%S')] $1${NC}"
}

error() {
    echo -e "${RED}[ERROR] $1${NC}"
}

warn() {
    echo -e "${YELLOW}[WARN] $1${NC}"
}

info() {
    echo -e "${BLUE}[INFO] $1${NC}"
}

# Check if port is available
check_port() {
    local port=$1
    if lsof -Pi :$port -sTCP:LISTEN -t >/dev/null 2>&1; then
        return 0
    else
        return 1
    fi
}

# Wait for port to be available
wait_for_port() {
    local port=$1
    local timeout=${2:-30}
    local count=0
    
    while ! check_port $port && [ $count -lt $timeout ]; do
        sleep 1
        ((count++))
    done
    
    if [ $count -eq $timeout ]; then
        error "Timeout waiting for port $port"
        return 1
    fi
    return 0
}

# Start nodemon with debugging
start_nodemon() {
    log "🚀 Starting nodemon with debugging..."
    
    # Check if nodemon is installed
    if ! command -v nodemon &> /dev/null; then
        error "nodemon not found. Installing..."
        npm install -g nodemon
    fi
    
    # Set environment variables
    export PY_DEBUG=1
    export NODE_ENV=development
    
    # Start nodemon with inspect
    nodemon --inspect=$NODE_DEBUG_PORT server.js &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES nodemon"
    
    # Wait for debugger to be ready
    log "⏳ Waiting for Node.js debugger on port $NODE_DEBUG_PORT..."
    wait_for_port $NODE_DEBUG_PORT
    
    log "✅ Nodemon started with debugging on port $NODE_DEBUG_PORT"
}

# Start Node.js server directly
start_node_server() {
    log "🚀 Starting Node.js server with debugging..."
    
    export PY_DEBUG=1
    export NODE_ENV=development
    
    node --inspect=$NODE_DEBUG_PORT server.js &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES node"
    
    wait_for_port $NODE_DEBUG_PORT
    log "✅ Node.js server started with debugging on port $NODE_DEBUG_PORT"
}

# Start Chrome with debugging
start_chrome_debug() {
    log "🌐 Starting Chrome with debugging..."
    
    # Check if Chrome is available
    local chrome_path=""
    if command -v google-chrome &> /dev/null; then
        chrome_path="google-chrome"
    elif command -v chromium-browser &> /dev/null; then
        chrome_path="chromium-browser"
    elif command -v /Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome &> /dev/null; then
        chrome_path="/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
    else
        error "Chrome not found. Please install Chrome or Chromium."
        return 1
    fi
    
    # Wait for server to be ready
    wait_for_port $NODE_PORT
    
    # Start Chrome with debugging flags
    $chrome_path \
        --disable-web-security \
        --disable-features=VizDisplayCompositor \
        --remote-debugging-port=$CHROME_DEBUG_PORT \
        --user-data-dir=/tmp/chrome-debug \
        http://localhost:$NODE_PORT &
    
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES chrome"
    
    wait_for_port $CHROME_DEBUG_PORT
    log "✅ Chrome started with debugging on port $CHROME_DEBUG_PORT"
    log "🌐 Application available at: http://localhost:$NODE_PORT"
}

# Start Edge with debugging
start_edge_debug() {
    log "🌐 Starting Edge with debugging..."
    
    local edge_path=""
    if command -v msedge &> /dev/null; then
        edge_path="msedge"
    elif command -v /Applications/Microsoft\ Edge.app/Contents/MacOS/Microsoft\ Edge &> /dev/null; then
        edge_path="/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge"
    else
        error "Edge not found. Please install Microsoft Edge."
        return 1
    fi
    
    wait_for_port $NODE_PORT
    
    $edge_path \
        --remote-debugging-port=$CHROME_DEBUG_PORT \
        --user-data-dir=/tmp/edge-debug \
        http://localhost:$NODE_PORT &
    
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES edge"
    
    wait_for_port $CHROME_DEBUG_PORT
    log "✅ Edge started with debugging on port $CHROME_DEBUG_PORT"
}

# Start Firefox with debugging
start_firefox_debug() {
    log "🦊 Starting Firefox with debugging..."
    
    if ! command -v firefox &> /dev/null; then
        error "Firefox not found. Please install Firefox."
        return 1
    fi
    
    wait_for_port $NODE_PORT
    
    firefox \
        --remote-debugging-port=$CHROME_DEBUG_PORT \
        --profile /tmp/firefox-debug \
        http://localhost:$NODE_PORT &
    
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES firefox"
    
    log "✅ Firefox started with debugging"
}

# Start Python debugging
start_python_debug() {
    log "🐍 Starting Python debugging..."
    
    # Check if debugpy is available
    if ! python -c "import debugpy" 2>/dev/null; then
        error "debugpy not found. Installing..."
        pip install debugpy
    fi
    
    export PYTHONPATH="$WORKSPACE_ROOT"
    
    # Start Python debug server
    python -m debugpy --listen 0.0.0.0:$PYTHON_DEBUG_PORT --wait-for-client &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES python_debug"
    
    wait_for_port $PYTHON_DEBUG_PORT
    log "✅ Python debug server started on port $PYTHON_DEBUG_PORT"
    log "🔗 Connect your Python debugger to localhost:$PYTHON_DEBUG_PORT"
}

# Run Python tests with debugging
run_python_tests() {
    log "🐍 Running Python tests with debugging..."
    
    export PYTHONPATH="$WORKSPACE_ROOT"
    
    python -m debugpy --listen 0.0.0.0:$PYTHON_DEBUG_PORT --wait-for-client -m pytest tests/test_convert.py -v &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES python_tests"
    
    wait_for_port $PYTHON_DEBUG_PORT
    log "✅ Python tests started with debugging on port $PYTHON_DEBUG_PORT"
}

# Run Python script with debugging
run_python_script() {
    local script=${1:-"sdf.py"}
    local args=${2:-"--help"}
    
    log "🐍 Running Python script: $script with debugging..."
    
    export PYTHONPATH="$WORKSPACE_ROOT"
    
    python -m debugpy --listen 0.0.0.0:$PYTHON_DEBUG_PORT --wait-for-client "$script" $args &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES python_script"
    
    wait_for_port $PYTHON_DEBUG_PORT
    log "✅ Python script started with debugging on port $PYTHON_DEBUG_PORT"
}

# Run Jest tests
run_jest_tests() {
    log "🧪 Running Jest tests..."
    
    export NODE_ENV=test
    
    npx jest --runInBand --detectOpenHandles --watchAll &
    PIDS="$PIDS $!"
    PROCESS_NAMES="$PROCESS_NAMES jest"
    
    log "✅ Jest tests started in watch mode"
}

# Attach to existing processes
attach_to_nodemon() {
    log "🔗 Attaching to existing nodemon process..."
    
    if ! check_port $NODE_DEBUG_PORT; then
        error "No nodemon process found on port $NODE_DEBUG_PORT"
        return 1
    fi
    
    log "✅ Ready to attach debugger to nodemon on port $NODE_DEBUG_PORT"
}

attach_to_chrome() {
    log "🔗 Attaching to existing Chrome process..."
    
    if ! check_port $CHROME_DEBUG_PORT; then
        error "No Chrome process found on port $CHROME_DEBUG_PORT"
        return 1
    fi
    
    log "✅ Ready to attach debugger to Chrome on port $CHROME_DEBUG_PORT"
}

attach_to_python() {
    log "🔗 Attaching to existing Python process..."
    
    if ! check_port $PYTHON_DEBUG_PORT; then
        error "No Python debug process found on port $PYTHON_DEBUG_PORT"
        return 1
    fi
    
    log "✅ Ready to attach debugger to Python on port $PYTHON_DEBUG_PORT"
}

# Show status
show_status() {
    echo -e "${CYAN}📊 Current Debug Status:${NC}"
    echo "=================================="
    
    if [ -n "$PIDS" ]; then
        local i=1
        for pid in $PIDS; do
            local name=$(echo $PROCESS_NAMES | cut -d' ' -f$i)
            if kill -0 "$pid" 2>/dev/null; then
                echo -e "${GREEN}✅ $name (PID: $pid)${NC}"
            else
                echo -e "${RED}❌ $name (PID: $pid) - Not running${NC}"
            fi
            i=$((i+1))
        done
    fi
    
    echo ""
    echo -e "${CYAN}🔌 Debug Ports:${NC}"
    echo "Node.js Debug: localhost:$NODE_DEBUG_PORT"
    echo "Chrome Debug: localhost:$CHROME_DEBUG_PORT"
    echo "Python Debug: localhost:$PYTHON_DEBUG_PORT"
    echo "App Server: localhost:$NODE_PORT"
}

# Main menu
show_menu() {
    echo -e "${PURPLE}🚀 Debug Launcher - VS Code Equivalent${NC}"
    echo "=============================================="
    echo ""
    echo -e "${CYAN}Compound Configurations:${NC}"
    echo "1) 🚀 Full Stack Hot Reload Debug (Nodemon + Chrome + Python)"
    echo "2) 🔗 Full Stack Attach Debug (Attach to existing processes)"
    echo "3) ⚡ Quick Debug (Backend + Frontend)"
    echo "4) 🐍 Python + Node Debug"
    echo ""
    echo -e "${CYAN}Individual Components:${NC}"
    echo "5) 🚀 Start Node.js Server with debugging"
    echo "6) 🔗 Attach to Nodemon"
    echo "7) 🌐 Start Chrome with debugging"
    echo "8) 🌐 Start Edge with debugging"
    echo "9) 🦊 Start Firefox with debugging"
    echo "10) 🔗 Attach to Chrome"
    echo "11) 🐍 Start Python debugging"
    echo "12) 🐍 Run Python tests with debugging"
    echo "13) 🐍 Run Python script with debugging"
    echo "14) 🔗 Attach to Python"
    echo "15) 🧪 Run Jest tests"
    echo ""
    echo -e "${CYAN}Utilities:${NC}"
    echo "16) 📊 Show status"
    echo "17) 🛑 Stop all processes"
    echo "18) ❌ Exit"
    echo ""
}

# Main function
main() {
    case "${1:-menu}" in
        "menu"|"")
            while true; do
                show_menu
                read -p "Select option (1-18): " choice
                
                case $choice in
                    1) full_stack_hot_reload ;;
                    2) full_stack_attach ;;
                    3) quick_debug ;;
                    4) python_node_debug ;;
                    5) start_node_server ;;
                    6) attach_to_nodemon ;;
                    7) start_chrome_debug ;;
                    8) start_edge_debug ;;
                    9) start_firefox_debug ;;
                    10) attach_to_chrome ;;
                    11) start_python_debug ;;
                    12) run_python_tests ;;
                    13) run_python_script ;;
                    14) attach_to_python ;;
                    15) run_jest_tests ;;
                    16) show_status ;;
                    17) cleanup ;;
                    18) exit 0 ;;
                    *) echo -e "${RED}Invalid option${NC}" ;;
                esac
                
                echo ""
                read -p "Press Enter to continue..."
                clear
            done
            ;;
        "full-stack-hot-reload"|"1")
            full_stack_hot_reload
            ;;
        "full-stack-attach"|"2")
            full_stack_attach
            ;;
        "quick-debug"|"3")
            quick_debug
            ;;
        "python-node"|"4")
            python_node_debug
            ;;
        "node-server"|"5")
            start_node_server
            ;;
        "attach-nodemon"|"6")
            attach_to_nodemon
            ;;
        "chrome"|"7")
            start_chrome_debug
            ;;
        "edge"|"8")
            start_edge_debug
            ;;
        "firefox"|"9")
            start_firefox_debug
            ;;
        "attach-chrome"|"10")
            attach_to_chrome
            ;;
        "python-debug"|"11")
            start_python_debug
            ;;
        "python-tests"|"12")
            run_python_tests
            ;;
        "python-script"|"13")
            run_python_script "${2:-sdf.py}" "${3:---help}"
            ;;
        "attach-python"|"14")
            attach_to_python
            ;;
        "jest"|"15")
            run_jest_tests
            ;;
        "status"|"16")
            show_status
            ;;
        "stop"|"17")
            cleanup
            ;;
        "help"|"-h"|"--help")
            echo "Usage: $0 [mode] [options]"
            echo ""
            echo "Modes:"
            echo "  menu                    Interactive menu (default)"
            echo "  full-stack-hot-reload   Start full stack with hot reload"
            echo "  full-stack-attach       Attach to existing processes"
            echo "  quick-debug             Quick backend + frontend debug"
            echo "  python-node             Python + Node debug"
            echo "  node-server             Start Node.js server"
            echo "  attach-nodemon          Attach to nodemon"
            echo "  chrome                  Start Chrome debug"
            echo "  edge                    Start Edge debug"
            echo "  firefox                 Start Firefox debug"
            echo "  attach-chrome           Attach to Chrome"
            echo "  python-debug            Start Python debug"
            echo "  python-tests            Run Python tests"
            echo "  python-script [script]  Run Python script"
            echo "  attach-python           Attach to Python"
            echo "  jest                    Run Jest tests"
            echo "  status                  Show status"
            echo "  stop                    Stop all processes"
            ;;
        *)
            error "Unknown mode: $1"
            echo "Use '$0 help' for usage information"
            exit 1
            ;;
    esac
}

# Compound configurations
full_stack_hot_reload() {
    log "🚀 Starting Full Stack Hot Reload Debug..."
    start_nodemon
    sleep 2
    start_chrome_debug
    sleep 1
    start_python_debug
    show_status
}

full_stack_attach() {
    log "🔗 Full Stack Attach Debug..."
    attach_to_nodemon
    attach_to_chrome
    attach_to_python
    show_status
}

quick_debug() {
    log "⚡ Quick Debug (Backend + Frontend)..."
    start_nodemon
    sleep 2
    start_chrome_debug
    show_status
}

python_node_debug() {
    log "🐍 Python + Node Debug..."
    start_nodemon
    sleep 2
    start_python_debug
    show_status
}

# Run main function
main "$@" 