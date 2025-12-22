#!/bin/bash
# Setup: prompts for API keys, saves to .env, installs dependencies
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/.env"

echo "=== pro_eval Setup ==="

# Prompt for API keys (hidden input)
read -r -s -p "OPENAI_API_KEY: " OPENAI_API_KEY && echo
read -r -s -p "TAMARIND_API_KEY: " TAMARIND_API_KEY && echo

# Validate keys are not empty
if [[ -z "$OPENAI_API_KEY" || -z "$TAMARIND_API_KEY" ]]; then
    echo "Error: API keys cannot be empty" >&2
    exit 1
fi

# Save to .env
cat > "$ENV_FILE" << EOF
OPENAI_API_KEY=$OPENAI_API_KEY
TAMARIND_API_KEY=$TAMARIND_API_KEY
EOF
echo "✓ Saved $ENV_FILE"

# Install dependencies
cd "$SCRIPT_DIR/tasks"
pip install -e .

echo "✓ Setup complete"
