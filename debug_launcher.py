#!/usr/bin/env python3
"""
Debug Launcher Script - Functional equivalent to VS Code launch.json & tasks.json
Usage: python debug_launcher.py [mode] [options]
"""

import os
import sys
import time
import signal
import subprocess
import threading
import socket
import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import argparse

# Colors for output
class Colors:
    RED = '\033[0;31m'
    GREEN = '\033[0;32m'
    YELLOW = '\033[1;33m'
    BLUE = '\033[0;34m'
    PURPLE = '\033[0;35m'
    CYAN = '\033[0;36m'
    NC = '\033[0m'  # No Color

class DebugLauncher:
    def __init__(self):
        """Initialize debug launcher with project paths and configuration"""
        self.workspace_root = os.getcwd()
        self.node_port = int(os.environ.get('PORT', 8080))  # Cloud Functions standard, configurable
        self.node_debug_port = 9229
        self.chrome_debug_port = 9223
        self.python_debug_port = 5678
        
        # Process management
        self.processes: Dict[str, subprocess.Popen] = {}
        self.process_names: Dict[str, str] = {}
        
        # Setup signal handlers
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
    
    def _signal_handler(self, signum, frame):
        """Handle cleanup on exit"""
        self.cleanup()
        sys.exit(0)
    
    def log(self, message: str):
        """Log with timestamp and color"""
        timestamp = time.strftime('%H:%M:%S')
        print(f"{Colors.GREEN}[{timestamp}] {message}{Colors.NC}")
    
    def error(self, message: str):
        """Error logging"""
        print(f"{Colors.RED}[ERROR] {message}{Colors.NC}")
    
    def warn(self, message: str):
        """Warning logging"""
        print(f"{Colors.YELLOW}[WARN] {message}{Colors.NC}")
    
    def info(self, message: str):
        """Info logging"""
        print(f"{Colors.BLUE}[INFO] {message}{Colors.NC}")
    
    def check_port(self, port: int) -> bool:
        """Check if port is available"""
        try:
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.settimeout(1)
                result = s.connect_ex(('localhost', port))
                return result == 0
        except:
            return False
    
    def wait_for_port(self, port: int, timeout: int = 30) -> bool:
        """Wait for port to be available"""
        start_time = time.time()
        while time.time() - start_time < timeout:
            if self.check_port(port):
                return True
            time.sleep(1)
        self.error(f"Timeout waiting for port {port}")
        return False
    
    def cleanup(self):
        """Clean up all running processes"""
        print(f"{Colors.YELLOW}🛑 Cleaning up processes...{Colors.NC}")
        for name, process in self.processes.items():
            if process.poll() is None:  # Process is still running
                print(f"{Colors.YELLOW}Killing process {name} (PID: {process.pid}){Colors.NC}")
                try:
                    process.terminate()
                    process.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    process.kill()
                except:
                    pass
    
    def start_nodemon(self):
        """Start nodemon with debugging"""
        self.log("🚀 Starting nodemon with debugging...")
        
        # Check if nodemon is installed
        try:
            subprocess.run(['nodemon', '--version'], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.error("nodemon not found. Installing...")
            subprocess.run(['npm', 'install', '-g', 'nodemon'], check=True)
        
        # Set environment variables
        env = os.environ.copy()
        env['PY_DEBUG'] = '1'
        env['NODE_ENV'] = 'development'
        
        # Start nodemon with inspect
        process = subprocess.Popen(
            ['nodemon', f'--inspect={self.node_debug_port}', 'server.js'],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['nodemon'] = process
        self.process_names['nodemon'] = 'Nodemon Server'
        
        # Wait for debugger to be ready
        self.log(f"⏳ Waiting for Node.js debugger on port {self.node_debug_port}...")
        if self.wait_for_port(self.node_debug_port):
            self.log(f"✅ Nodemon started with debugging on port {self.node_debug_port}")
        else:
            self.error("Failed to start nodemon")
    
    def start_node_server(self):
        """Start Node.js server directly"""
        self.log("🚀 Starting Node.js server with debugging...")
        
        env = os.environ.copy()
        env['PY_DEBUG'] = '1'
        env['NODE_ENV'] = 'development'
        
        process = subprocess.Popen(
            ['node', f'--inspect={self.node_debug_port}', 'server.js'],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['node'] = process
        self.process_names['node'] = 'Node.js Server'
        
        if self.wait_for_port(self.node_debug_port):
            self.log(f"✅ Node.js server started with debugging on port {self.node_debug_port}")
        else:
            self.error("Failed to start Node.js server")
    
    def start_chrome_debug(self):
        """Start Chrome with debugging"""
        self.log("🌐 Starting Chrome with debugging...")
        
        # Find Chrome executable
        chrome_paths = [
            'google-chrome',
            'chromium-browser',
            '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
            'C:\\Program Files\\Google\\Chrome\\Application\\chrome.exe',
            'C:\\Program Files (x86)\\Google\\Chrome\\Application\\chrome.exe'
        ]
        
        chrome_path = None
        for path in chrome_paths:
            try:
                subprocess.run([path, '--version'], check=True, capture_output=True)
                chrome_path = path
                break
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue
        
        if not chrome_path:
            self.error("Chrome not found. Please install Chrome or Chromium.")
            return
        
        # Wait for server to be ready
        self.wait_for_port(self.node_port)
        
        # Start Chrome with debugging flags
        chrome_args = [
            chrome_path,
            '--disable-web-security',
            '--disable-features=VizDisplayCompositor',
            f'--remote-debugging-port={self.chrome_debug_port}',
            '--user-data-dir=/tmp/chrome-debug',
            f'http://localhost:{self.node_port}'
        ]
        
        process = subprocess.Popen(
            chrome_args,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        
        self.processes['chrome'] = process
        self.process_names['chrome'] = 'Chrome Debug'
        
        if self.wait_for_port(self.chrome_debug_port):
            self.log(f"✅ Chrome started with debugging on port {self.chrome_debug_port}")
            self.log(f"🌐 Application available at: http://localhost:{self.node_port}")
        else:
            self.error("Failed to start Chrome")
    
    def start_edge_debug(self):
        """Start Edge with debugging"""
        self.log("🌐 Starting Edge with debugging...")
        
        edge_paths = [
            'msedge',
            '/Applications/Microsoft Edge.app/Contents/MacOS/Microsoft Edge',
            'C:\\Program Files (x86)\\Microsoft\\Edge\\Application\\msedge.exe'
        ]
        
        edge_path = None
        for path in edge_paths:
            try:
                subprocess.run([path, '--version'], check=True, capture_output=True)
                edge_path = path
                break
            except (subprocess.CalledProcessError, FileNotFoundError):
                continue
        
        if not edge_path:
            self.error("Edge not found. Please install Microsoft Edge.")
            return
        
        self.wait_for_port(self.node_port)
        
        edge_args = [
            edge_path,
            f'--remote-debugging-port={self.chrome_debug_port}',
            '--user-data-dir=/tmp/edge-debug',
            f'http://localhost:{self.node_port}'
        ]
        
        process = subprocess.Popen(
            edge_args,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        
        self.processes['edge'] = process
        self.process_names['edge'] = 'Edge Debug'
        
        if self.wait_for_port(self.chrome_debug_port):
            self.log(f"✅ Edge started with debugging on port {self.chrome_debug_port}")
        else:
            self.error("Failed to start Edge")
    
    def start_firefox_debug(self):
        """Start Firefox with debugging"""
        self.log("🦊 Starting Firefox with debugging...")
        
        try:
            subprocess.run(['firefox', '--version'], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            self.error("Firefox not found. Please install Firefox.")
            return
        
        self.wait_for_port(self.node_port)
        
        firefox_args = [
            'firefox',
            f'--remote-debugging-port={self.chrome_debug_port}',
            '--profile', '/tmp/firefox-debug',
            f'http://localhost:{self.node_port}'
        ]
        
        process = subprocess.Popen(
            firefox_args,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL
        )
        
        self.processes['firefox'] = process
        self.process_names['firefox'] = 'Firefox Debug'
        
        self.log("✅ Firefox started with debugging")
    
    def start_python_debug(self):
        """Start Python debugging"""
        self.log("🐍 Starting Python debugging...")
        
        # Check if debugpy is available
        try:
            import debugpy
        except ImportError:
            self.error("debugpy not found. Installing...")
            subprocess.run([sys.executable, '-m', 'pip', 'install', 'debugpy'], check=True)
        
        env = os.environ.copy()
        env['PYTHONPATH'] = str(self.workspace_root)
        
        process = subprocess.Popen(
            [sys.executable, '-m', 'debugpy', f'--listen=0.0.0.0:{self.python_debug_port}', '--wait-for-client'],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['python_debug'] = process
        self.process_names['python_debug'] = 'Python Debug Server'
        
        if self.wait_for_port(self.python_debug_port):
            self.log(f"✅ Python debug server started on port {self.python_debug_port}")
            self.log(f"🔗 Connect your Python debugger to localhost:{self.python_debug_port}")
        else:
            self.error("Failed to start Python debug server")
    
    def run_python_tests(self):
        """Run Python tests with debugging"""
        self.log("🐍 Running Python tests with debugging...")
        
        env = os.environ.copy()
        env['PYTHONPATH'] = str(self.workspace_root)
        
        process = subprocess.Popen(
            [sys.executable, '-m', 'debugpy', f'--listen=0.0.0.0:{self.python_debug_port}', '--wait-for-client', '-m', 'pytest', 'tests/test_convert.py', '-v'],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['python_tests'] = process
        self.process_names['python_tests'] = 'Python Tests'
        
        if self.wait_for_port(self.python_debug_port):
            self.log(f"✅ Python tests started with debugging on port {self.python_debug_port}")
        else:
            self.error("Failed to start Python tests")
    
    def run_python_script(self, script: str = "sdf.py", args: str = "--help"):
        """Run Python script with debugging"""
        self.log(f"🐍 Running Python script: {script} with debugging...")
        
        env = os.environ.copy()
        env['PYTHONPATH'] = str(self.workspace_root)
        
        script_args = args.split() if args else []
        cmd = [sys.executable, '-m', 'debugpy', f'--listen=0.0.0.0:{self.python_debug_port}', '--wait-for-client', script] + script_args
        
        process = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['python_script'] = process
        self.process_names['python_script'] = 'Python Script'
        
        if self.wait_for_port(self.python_debug_port):
            self.log(f"✅ Python script started with debugging on port {self.python_debug_port}")
        else:
            self.error("Failed to start Python script")
    
    def run_jest_tests(self):
        """Run Jest tests"""
        self.log("🧪 Running Jest tests...")
        
        env = os.environ.copy()
        env['NODE_ENV'] = 'test'
        
        process = subprocess.Popen(
            ['npx', 'jest', '--runInBand', '--detectOpenHandles', '--watchAll'],
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        self.processes['jest'] = process
        self.process_names['jest'] = 'Jest Tests'
        
        self.log("✅ Jest tests started in watch mode")
    
    def attach_to_nodemon(self):
        """Attach to existing nodemon process"""
        self.log("🔗 Attaching to existing nodemon process...")
        
        if not self.check_port(self.node_debug_port):
            self.error(f"No nodemon process found on port {self.node_debug_port}")
            return False
        
        self.log(f"✅ Ready to attach debugger to nodemon on port {self.node_debug_port}")
        return True
    
    def attach_to_chrome(self):
        """Attach to existing Chrome process"""
        self.log("🔗 Attaching to existing Chrome process...")
        
        if not self.check_port(self.chrome_debug_port):
            self.error(f"No Chrome process found on port {self.chrome_debug_port}")
            return False
        
        self.log(f"✅ Ready to attach debugger to Chrome on port {self.chrome_debug_port}")
        return True
    
    def attach_to_python(self):
        """Attach to existing Python process"""
        self.log("🔗 Attaching to existing Python process...")
        
        if not self.check_port(self.python_debug_port):
            self.error(f"No Python debug process found on port {self.python_debug_port}")
            return False
        
        self.log(f"✅ Ready to attach debugger to Python on port {self.python_debug_port}")
        return True
    
    def show_status(self):
        """Show current status"""
        print(f"{Colors.CYAN}📊 Current Debug Status:{Colors.NC}")
        print("==================================")
        
        for name, process in self.processes.items():
            process_name = self.process_names.get(name, name)
            if process.poll() is None:  # Process is still running
                print(f"{Colors.GREEN}✅ {process_name} (PID: {process.pid}){Colors.NC}")
            else:
                print(f"{Colors.RED}❌ {process_name} (PID: {process.pid}) - Not running{Colors.NC}")
        
        print("")
        print(f"{Colors.CYAN}🔌 Debug Ports:{Colors.NC}")
        print(f"Node.js Debug: localhost:{self.node_debug_port}")
        print(f"Chrome Debug: localhost:{self.chrome_debug_port}")
        print(f"Python Debug: localhost:{self.python_debug_port}")
        print(f"App Server: localhost:{self.node_port}")
    
    def show_menu(self):
        """Show interactive menu"""
        print(f"{Colors.PURPLE}🚀 Debug Launcher - VS Code Equivalent{Colors.NC}")
        print("==============================================")
        print("")
        print(f"{Colors.CYAN}Compound Configurations:{Colors.NC}")
        print("1) 🚀 Full Stack Hot Reload Debug (Nodemon + Chrome + Python)")
        print("2) 🔗 Full Stack Attach Debug (Attach to existing processes)")
        print("3) ⚡ Quick Debug (Backend + Frontend)")
        print("4) 🐍 Python + Node Debug")
        print("")
        print(f"{Colors.CYAN}Individual Components:{Colors.NC}")
        print("5) 🚀 Start Node.js Server with debugging")
        print("6) 🔗 Attach to Nodemon")
        print("7) 🌐 Start Chrome with debugging")
        print("8) 🌐 Start Edge with debugging")
        print("9) 🦊 Start Firefox with debugging")
        print("10) 🔗 Attach to Chrome")
        print("11) 🐍 Start Python debugging")
        print("12) 🐍 Run Python tests with debugging")
        print("13) 🐍 Run Python script with debugging")
        print("14) 🔗 Attach to Python")
        print("15) 🧪 Run Jest tests")
        print("")
        print(f"{Colors.CYAN}Utilities:{Colors.NC}")
        print("16) 📊 Show status")
        print("17) 🛑 Stop all processes")
        print("18) ❌ Exit")
        print("")
    
    def interactive_menu(self):
        """Run interactive menu"""
        while True:
            self.show_menu()
            try:
                choice = input("Select option (1-18): ").strip()
                
                if choice == '1':
                    self.full_stack_hot_reload()
                elif choice == '2':
                    self.full_stack_attach()
                elif choice == '3':
                    self.quick_debug()
                elif choice == '4':
                    self.python_node_debug()
                elif choice == '5':
                    self.start_node_server()
                elif choice == '6':
                    self.attach_to_nodemon()
                elif choice == '7':
                    self.start_chrome_debug()
                elif choice == '8':
                    self.start_edge_debug()
                elif choice == '9':
                    self.start_firefox_debug()
                elif choice == '10':
                    self.attach_to_chrome()
                elif choice == '11':
                    self.start_python_debug()
                elif choice == '12':
                    self.run_python_tests()
                elif choice == '13':
                    script = input("Enter script name (default: sdf.py): ").strip() or "sdf.py"
                    args = input("Enter arguments (default: --help): ").strip() or "--help"
                    self.run_python_script(script, args)
                elif choice == '14':
                    self.attach_to_python()
                elif choice == '15':
                    self.run_jest_tests()
                elif choice == '16':
                    self.show_status()
                elif choice == '17':
                    self.cleanup()
                elif choice == '18':
                    break
                else:
                    print(f"{Colors.RED}Invalid option{Colors.NC}")
                
                print("")
                input("Press Enter to continue...")
                os.system('clear' if os.name == 'posix' else 'cls')
                
            except KeyboardInterrupt:
                break
    
    # Compound configurations
    def full_stack_hot_reload(self):
        """Full stack hot reload debug"""
        self.log("🚀 Starting Full Stack Hot Reload Debug...")
        self.start_nodemon()
        time.sleep(2)
        self.start_chrome_debug()
        time.sleep(1)
        self.start_python_debug()
        self.show_status()
    
    def full_stack_attach(self):
        """Full stack attach debug"""
        self.log("🔗 Full Stack Attach Debug...")
        self.attach_to_nodemon()
        self.attach_to_chrome()
        self.attach_to_python()
        self.show_status()
    
    def quick_debug(self):
        """Quick debug (backend + frontend)"""
        self.log("⚡ Quick Debug (Backend + Frontend)...")
        self.start_nodemon()
        time.sleep(2)
        self.start_chrome_debug()
        self.show_status()
    
    def python_node_debug(self):
        """Python + Node debug"""
        self.log("🐍 Python + Node Debug...")
        self.start_nodemon()
        time.sleep(2)
        self.start_python_debug()
        self.show_status()

def main():
    parser = argparse.ArgumentParser(description='Debug Launcher - VS Code Equivalent')
    parser.add_argument('mode', nargs='?', default='menu', help='Debug mode to run')
    parser.add_argument('--script', default='sdf.py', help='Python script to run (for python-script mode)')
    parser.add_argument('--args', default='--help', help='Arguments for Python script')
    
    args = parser.parse_args()
    
    launcher = DebugLauncher()
    
    try:
        if args.mode == 'menu':
            launcher.interactive_menu()
        elif args.mode == 'full-stack-hot-reload':
            launcher.full_stack_hot_reload()
        elif args.mode == 'full-stack-attach':
            launcher.full_stack_attach()
        elif args.mode == 'quick-debug':
            launcher.quick_debug()
        elif args.mode == 'python-node':
            launcher.python_node_debug()
        elif args.mode == 'node-server':
            launcher.start_node_server()
        elif args.mode == 'attach-nodemon':
            launcher.attach_to_nodemon()
        elif args.mode == 'chrome':
            launcher.start_chrome_debug()
        elif args.mode == 'edge':
            launcher.start_edge_debug()
        elif args.mode == 'firefox':
            launcher.start_firefox_debug()
        elif args.mode == 'attach-chrome':
            launcher.attach_to_chrome()
        elif args.mode == 'python-debug':
            launcher.start_python_debug()
        elif args.mode == 'python-tests':
            launcher.run_python_tests()
        elif args.mode == 'python-script':
            launcher.run_python_script(args.script, args.args)
        elif args.mode == 'attach-python':
            launcher.attach_to_python()
        elif args.mode == 'jest':
            launcher.run_jest_tests()
        elif args.mode == 'status':
            launcher.show_status()
        elif args.mode == 'stop':
            launcher.cleanup()
        else:
            launcher.error(f"Unknown mode: {args.mode}")
            print("Available modes:")
            print("  menu                    Interactive menu (default)")
            print("  full-stack-hot-reload   Start full stack with hot reload")
            print("  full-stack-attach       Attach to existing processes")
            print("  quick-debug             Quick backend + frontend debug")
            print("  python-node             Python + Node debug")
            print("  node-server             Start Node.js server")
            print("  attach-nodemon          Attach to nodemon")
            print("  chrome                  Start Chrome debug")
            print("  edge                    Start Edge debug")
            print("  firefox                 Start Firefox debug")
            print("  attach-chrome           Attach to Chrome")
            print("  python-debug            Start Python debug")
            print("  python-tests            Run Python tests")
            print("  python-script           Run Python script")
            print("  attach-python           Attach to Python")
            print("  jest                    Run Jest tests")
            print("  status                  Show status")
            print("  stop                    Stop all processes")
            sys.exit(1)
    
    except KeyboardInterrupt:
        print("\nInterrupted by user")
    finally:
        launcher.cleanup()

if __name__ == '__main__':
    main() 