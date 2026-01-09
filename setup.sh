#!/bin/bash
# Setup: checks env vars for API keys, saves to .env, installs dependencies
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="$SCRIPT_DIR/.env"

echo "=== pro_eval Setup ==="

# Validate keys are present in environment
if [[ -z "$OPENAI_API_KEY" ]]; then
    echo "Error: OPENAI_API_KEY environment variable is not set" >&2
    exit 1
fi

if [[ -z "$TAMARIND_API_KEY" ]]; then
    echo "Error: TAMARIND_API_KEY environment variable is not set" >&2
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
