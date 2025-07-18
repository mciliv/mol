#!/bin/bash

echo "Setting up mol project..."

# Install dependencies
npm install

# Setup aliases
./setup-aliases.sh setup

echo "Setup complete!"
echo "Run 'source ~/.zshrc' (or ~/.bashrc) to load aliases"
echo "Then use commands like: dev, test, deploy, etc." 