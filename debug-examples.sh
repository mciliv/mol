#!/bin/bash

# Debug Launcher Examples
# This script shows examples of how to use the debug launcher

echo "🚀 Debug Launcher Examples"
echo "=========================="
echo ""

echo "📋 Available Commands:"
echo ""

echo "🔧 Interactive Mode:"
echo "  ./debug                    # Interactive menu (bash)"
echo "  ./debug --python          # Interactive menu (Python)"
echo ""

echo "🚀 Compound Configurations (VS Code Compounds):"
echo "  ./debug full-stack-hot-reload    # Nodemon + Chrome + Python"
echo "  ./debug full-stack-attach        # Attach to existing processes"
echo "  ./debug quick-debug              # Nodemon + Chrome"
echo "  ./debug python-node              # Nodemon + Python"
echo ""

echo "🔧 Individual Components:"
echo "  ./debug node-server               # Start Node.js server"
echo "  ./debug attach-nodemon            # Attach to nodemon"
echo "  ./debug chrome                    # Start Chrome debug"
echo "  ./debug edge                      # Start Edge debug"
echo "  ./debug firefox                   # Start Firefox debug"
echo "  ./debug attach-chrome             # Attach to Chrome"
echo "  ./debug python-debug              # Start Python debug"
echo "  ./debug python-tests              # Run Python tests"
echo "  ./debug python-script sdf.py      # Run Python script"
echo "  ./debug attach-python             # Attach to Python"
echo "  ./debug jest                      # Run Jest tests"
echo ""

echo "🛠️  Utilities:"
echo "  ./debug status                    # Show status"
echo "  ./debug stop                      # Stop all processes"
echo "  ./debug help                      # Show help"
echo ""

echo "🐍 Python Version Examples:"
echo "  ./debug --python full-stack-hot-reload"
echo "  ./debug --python python-script --script uniprot.py --args '--version'"
echo ""

echo "🔌 Debug Ports:"
echo "  Node.js Debug: localhost:9229"
echo "  Chrome Debug: localhost:9223"
echo "  Python Debug: localhost:5678"
echo "  App Server: localhost:3000"
echo ""

echo "📖 For detailed documentation, see: DEBUG_LAUNCHER.md"
echo ""

echo "💡 Quick Start Examples:"
echo ""

echo "1. Start full stack development:"
echo "   ./debug full-stack-hot-reload"
echo ""

echo "2. Quick backend + frontend debug:"
echo "   ./debug quick-debug"
echo ""

echo "3. Python debugging only:"
echo "   ./debug python-debug"
echo ""

echo "4. Run Python tests with debugging:"
echo "   ./debug python-tests"
echo ""

echo "5. Check what's running:"
echo "   ./debug status"
echo ""

echo "6. Stop all debug processes:"
echo "   ./debug stop"
echo "" 