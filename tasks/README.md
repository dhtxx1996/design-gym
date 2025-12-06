# Computational Biology Tasks

AI agent for computational biology tasks using Tamarind Bio API.

## Quick Setup

```bash
# Install
pip install -e .

# Set API keys
cp env.example .env
# Edit .env with your keys:
#   TAMARIND_API_KEY=your-key-from-tamarind.bio
#   OPENAI_API_KEY=your-openai-key
```

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
