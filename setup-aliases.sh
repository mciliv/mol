#!/bin/bash

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RC_FILE=~/.zshrc
[[ -f ~/.bashrc ]] && RC_FILE=~/.bashrc

# Remove existing aliases and add new ones
sed -i.bak '/# mol aliases/,/# end mol aliases/d' "$RC_FILE"

cat >> "$RC_FILE" << EOF

# mol aliases
alias dev="$PROJECT_DIR/run dev"
alias test="$PROJECT_DIR/run test"
alias deploy="$PROJECT_DIR/run deploy"
alias format="$PROJECT_DIR/run format"
alias server="$PROJECT_DIR/run server"
alias debug="$PROJECT_DIR/run debug"
alias unit="$PROJECT_DIR/run unit"
alias integration="$PROJECT_DIR/run integration"
alias system="$PROJECT_DIR/run system"
alias watch="$PROJECT_DIR/run watch"
alias pytest="$PROJECT_DIR/run pytest"
alias ship="$PROJECT_DIR/run ship"
alias ip="$PROJECT_DIR/run ip"
alias mobile="$PROJECT_DIR/run mobile"
alias cert="$PROJECT_DIR/run cert"
alias commit="$PROJECT_DIR/run commit"
# end mol aliases
EOF

echo "Aliases added to $RC_FILE"
echo "Run 'source $RC_FILE' to reload" 