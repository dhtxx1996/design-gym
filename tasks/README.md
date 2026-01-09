# Computational Biology Tasks

AI agent for computational biology tasks using Tamarind Bio API.

## Quick Setup

```bash
# Install
pip install -e .

# Set API keys as environment variables
# On Unix/Linux/Mac:
export TAMARIND_API_KEY=your-key-from-tamarind.bio
export OPENAI_API_KEY=your-openai-key

# On Windows (PowerShell):
# $env:TAMARIND_API_KEY="your-key-from-tamarind.bio"
# $env:OPENAI_API_KEY="your-openai-key"

# On Windows (CMD):
# set TAMARIND_API_KEY=your-key-from-tamarind.bio
# set OPENAI_API_KEY=your-openai-key
```

**Note:** Credentials are loaded from environment variables only (no `.env` file needed). This makes cross-platform setup easier. For persistent setup, add these to your shell profile or system environment variables.

## Run Agent

```bash
python agent.py --task ph_sensitive_design
```

## Create New Tasks

1. Create a folder: `my_task/`
2. Add `question.md` with task description
3. Add `data/` folder with input files
4. Run: `python agent.py --task my_task`

## Tamarind CLI (Optional)

```bash
tamarind --list-tools          # List available tools
tamarind --tool-info esmfold   # Get tool parameters
tamarind --test-esmfold        # Test with a sample sequence
```
