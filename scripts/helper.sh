#!/bin/sh

set -e

# Load configuration
source "$(dirname "$0")/config.sh"

PYVERSION=$PYTHON_VERSION

if ! command -v pyenv &> /dev/null; then
    curl https://pyenv.run | bash
    export PATH="$HOME/.pyenv/bin:$PATH"
    eval "$(pyenv init --path)"
    eval "$(pyenv init -)"
    eval "$(pyenv virtualenv-init -)"
fi
pyenv install -s $PYVERSION
pyenv local $PYVERSION

curl -sSL https://install.python-poetry.org | python3 -
export PATH="$HOME/.local/bin:$PATH"
poetry env use $PYVERSION
eval $(poetry env activate)
export PYTHONPATH=$(pwd)
poetry sync --no-root
echo "To deactivate env, run: exit"