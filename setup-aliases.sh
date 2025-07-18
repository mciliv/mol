#!/bin/bash

# Get the project directory
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create aliases for all run commands
setup_aliases() {
  echo "Setting up aliases for mol project..."
  
  # Check if .bashrc or .zshrc exists
  if [[ -f ~/.bashrc ]]; then
    RC_FILE=~/.bashrc
  elif [[ -f ~/.zshrc ]]; then
    RC_FILE=~/.zshrc
  else
    echo "No .bashrc or .zshrc found. Please manually add aliases to your shell config."
    return 1
  fi
  
  # Remove existing mol aliases
  sed -i.bak '/# mol aliases/,/# end mol aliases/d' "$RC_FILE"
  
  # Add new aliases
  cat >> "$RC_FILE" << EOF

# mol aliases
alias dev="$PROJECT_DIR/run dev"
alias start="$PROJECT_DIR/run start"
alias server="$PROJECT_DIR/run server"
alias debug="$PROJECT_DIR/run debug"
alias unsafe="$PROJECT_DIR/run unsafe"
alias test="$PROJECT_DIR/run test"
alias unit="$PROJECT_DIR/run unit"
alias integration="$PROJECT_DIR/run integration"
alias system="$PROJECT_DIR/run system"
alias watch="$PROJECT_DIR/run watch"
alias pytest="$PROJECT_DIR/run pytest"
alias deploy="$PROJECT_DIR/run deploy"
alias ship="$PROJECT_DIR/run ship"
alias format="$PROJECT_DIR/run format"
alias ip="$PROJECT_DIR/run ip"
alias mobile="$PROJECT_DIR/run mobile"
alias cert="$PROJECT_DIR/run cert"
alias tunnel="$PROJECT_DIR/run tunnel"
alias build="$PROJECT_DIR/run build"
alias commit="$PROJECT_DIR/run commit"
alias check-commit="$PROJECT_DIR/run check-commit"
alias domain-status="$PROJECT_DIR/run domain:status"
alias domain-setup="$PROJECT_DIR/run domain:setup"
# end mol aliases
EOF

  echo "Aliases added to $RC_FILE"
  echo "Run 'source $RC_FILE' to reload your shell config"
}

# Update aliases (same as setup)
update_aliases() {
  setup_aliases
}

case "$1" in
  "setup"|"")
    setup_aliases
    ;;
  "update")
    update_aliases
    ;;
  "help"|"-h"|"--help")
    echo "Usage: $0 [setup|update|help]"
    echo ""
    echo "  setup   - Set up aliases in shell config"
    echo "  update  - Update existing aliases"
    echo "  help    - Show this help"
    ;;
  *)
    echo "Unknown command: $1"
    echo "Run '$0 help' for available commands"
    exit 1
    ;;
esac 