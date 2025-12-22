# design-gym

A framework for evaluating AI agents on computational protein engineering tasks.

## Overview

This project provides:
- **Task-agnostic AI agent** that solves computational biology problems using tools like ProteinMPNN, ESMFold, and AlphaFold
- **Automated evaluation** system that scores agent outputs against rubrics
- **Benchmark tasks** derived from research papers, covering protein design, analysis, and engineering challenges

## Quick Start

```bash
# 1. Setup
cd tasks/
pip install -e .

# 2. Set API keys as environment variables
# On Unix/Linux/Mac:
export TAMARIND_API_KEY=your_tamarind_api_key_here
export OPENAI_API_KEY=your_openai_api_key_here

# On Windows (PowerShell):
# $env:TAMARIND_API_KEY="your_tamarind_api_key_here"
# $env:OPENAI_API_KEY="your_openai_api_key_here"

# On Windows (CMD):
# set TAMARIND_API_KEY=your_tamarind_api_key_here
# set OPENAI_API_KEY=your_openai_api_key_here

# 3. Run agent
python agent.py --task ph_sensitive_design --question question_easy.md --model gpt-4o --output my_run --overwrite
# Or: ./bench.sh easy gpt-4o my_run

# 4. Evaluate
python eval.py --task ph_sensitive_design --outputs my_run --k 3
```

## Structure

- `tasks/` - Agent code (`agent.py`, `eval.py`) and task directories
- `context/` - Question creation guides and examples
- `papers/` - Source research papers

## Features

- Tool integration (Tamarind Bio: ProteinMPNN, ESMFold, AlphaFold)
- Execution DAG tracking (visualize tool calls and data flow)
- LLM-based evaluation with rubric scoring
- Automatic result visualization

## Documentation

- Agent usage: `tasks/README.md`
- Evaluation workflow: `tasks/RUN_AND_EVAL.md`
- Creating tasks: `context/paperbench_guide.md`

## Output Files

After running and evaluating:
- `tasks/{task}/eval_results.json` - Evaluation scores
- `tasks/{task}/eval_results.png` - Score visualizations
- `tasks/{task}/agent_output/{run}/execution_dag.png` - Execution flow graph

## Requirements

- Python 3.8+
- OpenAI API key (for agent and evaluator) - set as `OPENAI_API_KEY` environment variable
- Tamarind Bio API key (for protein design tools) - set as `TAMARIND_API_KEY` environment variable
- See `tasks/requirements.txt` for Python dependencies

### Setting Environment Variables

Credentials are loaded from environment variables only (no `.env` file needed). This makes cross-platform setup easier:

**Unix/Linux/Mac (bash/zsh):**
```bash
export TAMARIND_API_KEY=your_key_here
export OPENAI_API_KEY=your_key_here
```

For persistent setup, add these to your shell profile (`.bashrc`, `.zshrc`, etc.) or system environment variables.
