#!/usr/bin/env python3
"""
Task-Agnostic AI Agent for Computational Biology

This agent uses OpenAI's function calling to work on computational biology
tasks defined in task directories. Each task should have:
- question.md: Task description and requirements
- data/: Input data files

Usage:
    python agent.py --task <task_name> [--max-iterations 20]
"""

import os
import sys
import json
import argparse
import traceback
from pathlib import Path
from datetime import datetime
from typing import Optional

from dotenv import load_dotenv
from openai import OpenAI

# Load environment variables
load_dotenv()
load_dotenv(Path(__file__).parent / ".env")

# Add current directory for imports
sys.path.insert(0, str(Path(__file__).parent))
from tamarind_client import TamarindClient


# =============================================================================
# Tool Definitions for OpenAI Function Calling
# =============================================================================

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "read_file",
            "description": "Read the contents of a file (PDB, JSON, CSV, Python, etc.).",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Path to the file (relative to task directory)"
                    }
                },
                "required": ["path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "write_file",
            "description": "Write content to a file in the output directory.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Path to write (relative to output directory)"
                    },
                    "content": {
                        "type": "string",
                        "description": "Content to write"
                    }
                },
                "required": ["path", "content"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "list_directory",
            "description": "List files and directories in a path.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Directory path (relative to task directory)"
                    }
                },
                "required": ["path"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "run_python",
            "description": "Execute Python code. Available: numpy (np), pandas (pd), BioPython (Bio, PDBParser, ShrakeRupley, seq1, SeqIO), json, Path, os. Use print() for output.",
            "parameters": {
                "type": "object",
                "properties": {
                    "code": {
                        "type": "string",
                        "description": "Python code to execute"
                    }
                },
                "required": ["code"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "tamarind_list_tools",
            "description": "List available tools on Tamarind Bio API (proteinmpnn, esmfold, alphafold, etc.).",
            "parameters": {"type": "object", "properties": {}, "required": []}
        }
    },
    {
        "type": "function",
        "function": {
            "name": "tamarind_get_tool_spec",
            "description": "Get the specification and parameters for a Tamarind tool. ALWAYS call this before using a tool to understand its parameters.",
            "parameters": {
                "type": "object",
                "properties": {
                    "tool_name": {
                        "type": "string",
                        "description": "Name of the tool (e.g., 'proteinmpnn', 'esmfold')"
                    }
                },
                "required": ["tool_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "tamarind_upload_file",
            "description": "Upload a file to Tamarind Bio (required before using it in jobs).",
            "parameters": {
                "type": "object",
                "properties": {
                    "filepath": {
                        "type": "string",
                        "description": "Path to file (relative to task directory)"
                    }
                },
                "required": ["filepath"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "tamarind_submit_job",
            "description": "Submit a job to Tamarind Bio and wait for results. Pass all parameters inside the 'params' dict.",
            "parameters": {
                "type": "object",
                "properties": {
                    "tool_name": {
                        "type": "string",
                        "description": "Tool name (e.g., 'proteinmpnn', 'esmfold')"
                    },
                    "params": {
                        "type": "object",
                        "description": "Job parameters as key-value pairs. Check tamarind_get_tool_spec for required params."
                    }
                },
                "required": ["tool_name", "params"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "task_complete",
            "description": "Call when you have completed the task with a summary.",
            "parameters": {
                "type": "object",
                "properties": {
                    "summary": {
                        "type": "string",
                        "description": "Summary of completed work and key findings"
                    }
                },
                "required": ["summary"]
            }
        }
    }
]


# =============================================================================
# Tool Implementations
# =============================================================================

class AgentTools:
    """Tool implementations for the agent with session-scoped job tracking."""
    
    def __init__(self, task_dir: Path, output_dir: Path, session_id: str):
        self.task_dir = task_dir
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.session_id = session_id
        self._tamarind: Optional[TamarindClient] = None
        
        # Track agent execution metadata
        self.execution_start = datetime.now()
        self.tool_calls: list[dict] = []
        
    @property
    def tamarind(self) -> TamarindClient:
        if self._tamarind is None:
            # Create client with session ID to track jobs
            self._tamarind = TamarindClient(session_id=self.session_id)
            print(f"[agent] Tamarind session: {self.session_id}")
        return self._tamarind
    
    def get_execution_summary(self) -> dict:
        """Get summary of agent execution for logging."""
        tamarind_summary = self.tamarind.get_session_summary() if self._tamarind else {}
        return {
            "session_id": self.session_id,
            "execution_start": self.execution_start.isoformat(),
            "execution_end": datetime.now().isoformat(),
            "tool_calls_count": len(self.tool_calls),
            "tamarind": tamarind_summary,
        }
    
    def read_file(self, path: str) -> str:
        full_path = self.task_dir / path
        if not full_path.exists():
            full_path = self.output_dir / path
        if not full_path.exists():
            return f"Error: File not found: {path}"
        
        try:
            content = full_path.read_text()
            if len(content) > 50000:
                return content[:50000] + f"\n\n... [truncated, {len(content)} total chars]"
            return content
        except Exception as e:
            return f"Error reading file: {e}"
    
    def write_file(self, path: str, content: str) -> str:
        full_path = self.output_dir / path
        full_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            full_path.write_text(content)
            return f"Successfully wrote {len(content)} bytes to {full_path}"
        except Exception as e:
            return f"Error writing file: {e}"
    
    def list_directory(self, path: str) -> str:
        full_path = self.task_dir / path
        if not full_path.exists():
            full_path = self.output_dir / path
        if not full_path.exists():
            return f"Error: Directory not found: {path}"
        
        try:
            items = []
            for item in sorted(full_path.iterdir()):
                prefix = "[DIR] " if item.is_dir() else "[FILE]"
                size = f"({item.stat().st_size} bytes)" if item.is_file() else ""
                items.append(f"{prefix} {item.name} {size}")
            return "\n".join(items) if items else "(empty directory)"
        except Exception as e:
            return f"Error listing directory: {e}"
    
    def run_python(self, code: str) -> str:
        import io
        from contextlib import redirect_stdout, redirect_stderr
        
        original_cwd = os.getcwd()
        os.chdir(self.task_dir)
        
        exec_globals = {
            "__builtins__": __builtins__,
            "json": json, "Path": Path, "os": os,
        }
        
        # Import libraries
        try:
            import numpy as np
            exec_globals["np"] = exec_globals["numpy"] = np
        except ImportError:
            pass
        
        try:
            import pandas as pd
            exec_globals["pd"] = exec_globals["pandas"] = pd
        except ImportError:
            pass
        
        try:
            import Bio
            from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Selection
            from Bio.PDB.SASA import ShrakeRupley
            from Bio.SeqUtils import seq1
            from Bio import SeqIO
            exec_globals.update({
                "Bio": Bio, "PDBParser": PDBParser, "ShrakeRupley": ShrakeRupley,
                "NeighborSearch": NeighborSearch, "PDBIO": PDBIO, "Selection": Selection,
                "seq1": seq1, "SeqIO": SeqIO
            })
        except ImportError:
            pass
        
        exec_globals["TASK_DIR"] = self.task_dir
        exec_globals["OUTPUT_DIR"] = self.output_dir
        exec_globals["task_dir"] = str(self.task_dir)
        exec_globals["output_dir"] = str(self.output_dir)
        
        stdout = io.StringIO()
        stderr = io.StringIO()
        
        try:
            with redirect_stdout(stdout), redirect_stderr(stderr):
                exec(code, exec_globals)
            
            output = stdout.getvalue()
            errors = stderr.getvalue()
            
            result = ""
            if output:
                result += f"Output:\n{output}"
            if errors:
                result += f"\nStderr:\n{errors}"
            if not result:
                result = "(Code executed successfully with no output)"
            
            if len(result) > 10000:
                result = result[:10000] + "\n... [truncated]"
            
            return result
            
        except Exception as e:
            return f"Error executing code:\n{traceback.format_exc()}"
        finally:
            os.chdir(original_cwd)
    
    def tamarind_list_tools(self) -> str:
        try:
            tools = self.tamarind.list_tool_names()
            return "Available tools:\n" + "\n".join(f"  - {t}" for t in tools)
        except Exception as e:
            return f"Error: {e}"
    
    def tamarind_get_tool_spec(self, tool_name: str) -> str:
        try:
            spec = self.tamarind.get_tool_spec(tool_name)
            if spec:
                return json.dumps(spec, indent=2)
            return f"Tool '{tool_name}' not found"
        except Exception as e:
            return f"Error: {e}"
    
    def tamarind_upload_file(self, filepath: str) -> str:
        full_path = self.task_dir / filepath
        if not full_path.exists():
            full_path = self.output_dir / filepath
        if not full_path.exists():
            return f"Error: File not found: {filepath}"
        
        try:
            result = self.tamarind.upload_file(str(full_path))
            # Provide clear instructions about how to reference the file
            filename = result.get("filename", Path(filepath).name)
            return (
                f"Upload successful!\n"
                f"  Filename on Tamarind: {filename}\n"
                f"  Session ID: {self.session_id}\n"
                f"  Files uploaded this session: {list(self.tamarind.get_session_files().keys())}\n\n"
                f"To use this file in a job, reference it as: \"{filename}\"\n"
                f"Example: {{\"pdbFile\": \"{filename}\"}}"
            )
        except Exception as e:
            return f"Error uploading: {e}"
    
    def tamarind_submit_job(self, tool_name: str, params: dict) -> str:
        try:
            result = self.tamarind.submit_job_sync(tool_name, params, timeout=600)
            
            if "job_name" in result:
                try:
                    download_path = self.tamarind.download_results(
                        result["job_name"], 
                        output_dir=str(self.output_dir / "tamarind_results")
                    )
                    result["downloaded_to"] = str(download_path)
                except Exception as download_error:
                    result["download_error"] = str(download_error)
                    result["note"] = "Job completed but result download failed. Results may be available on Tamarind web UI."
            
            return json.dumps(result, indent=2, default=str)
        except Exception as e:
            error_msg = str(e)
            
            # Provide helpful context based on error type
            if "500" in error_msg and "/jobs" in error_msg:
                return (
                    f"Error: Tamarind API status endpoint returned 500 error.\n"
                    f"This is a transient server issue on Tamarind's side.\n"
                    f"The job may have been submitted successfully - check https://app.tamarind.bio\n\n"
                    f"Details: {error_msg}"
                )
            elif "400" in error_msg:
                return (
                    f"Error: Bad request (400) - check your parameters.\n"
                    f"Common issues:\n"
                    f"  - File not uploaded (use tamarind_upload_file first)\n"
                    f"  - Wrong parameter format (dict vs string)\n"
                    f"  - Missing required parameters\n\n"
                    f"Details: {error_msg}"
                )
            else:
                return f"Error: {error_msg}"
    
    def execute_tool(self, tool_name: str, arguments: dict) -> str:
        handlers = {
            "read_file": lambda: self.read_file(arguments["path"]),
            "write_file": lambda: self.write_file(arguments["path"], arguments["content"]),
            "list_directory": lambda: self.list_directory(arguments["path"]),
            "run_python": lambda: self.run_python(arguments["code"]),
            "tamarind_list_tools": lambda: self.tamarind_list_tools(),
            "tamarind_get_tool_spec": lambda: self.tamarind_get_tool_spec(arguments["tool_name"]),
            "tamarind_upload_file": lambda: self.tamarind_upload_file(arguments["filepath"]),
            "tamarind_submit_job": lambda: self._handle_tamarind_submit(arguments),
            "task_complete": lambda: f"TASK_COMPLETE: {arguments['summary']}",
        }
        return handlers.get(tool_name, lambda: f"Unknown tool: {tool_name}")()
    
    def _handle_tamarind_submit(self, arguments: dict) -> str:
        tool = arguments.get("tool_name")
        params = dict(arguments.get("params", {}))
        # Merge extra top-level args into params (handles common mistake)
        for key in arguments:
            if key not in ["tool_name", "params"]:
                params[key] = arguments[key]
        return self.tamarind_submit_job(tool, params)


# =============================================================================
# Agent Loop
# =============================================================================

def run_agent(task_name: str, max_iterations: int = 20, model: str = "gpt-4o", output_name: str = None):
    """Run the agent on a task."""
    tasks_dir = Path(__file__).parent
    task_dir = tasks_dir / task_name
    name = output_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = task_dir / "agent_output" / name
    
    # Generate unique session ID for tracking jobs
    session_id = f"{name}_{datetime.now().strftime('%H%M%S_%f')}"
    
    if not task_dir.exists():
        print(f"Error: Task directory not found: {task_dir}")
        return
    
    # Ensure clean output directory (no contamination from previous runs)
    if output_dir.exists():
        print(f"Warning: Output directory already exists: {output_dir}")
        print("  Results from previous runs may contaminate evaluation.")
        print("  Consider using a fresh output name or clearing the directory.")
    
    print("=" * 60)
    print("AI Agent for Computational Biology")
    print("=" * 60)
    print(f"Task: {task_name}")
    print(f"Output: {output_dir}")
    print(f"Session ID: {session_id}")
    print(f"Model: {model}")
    print(f"Max iterations: {max_iterations}")
    print()
    
    client = OpenAI()
    tools = AgentTools(task_dir, output_dir, session_id)
    
    # Load task description from question.md if it exists
    question_path = task_dir / "question.md"
    if question_path.exists():
        task_description = question_path.read_text()
    else:
        task_description = "Complete the computational biology task in this directory. Explore the files to understand requirements."
    
    # Task-agnostic system prompt
    system_prompt = f"""You are an expert computational biologist. Complete the task step by step using the available tools.

TASK DESCRIPTION:
{task_description}

AVAILABLE TOOLS:
- read_file: Read files (PDB, JSON, CSV, Python, etc.)
- write_file: Save outputs to the output directory
- list_directory: Explore available files
- run_python: Execute Python code (numpy, pandas, BioPython available)
- tamarind_list_tools: List available Tamarind Bio ML tools
- tamarind_get_tool_spec: Get tool parameters (ALWAYS call before using a tool)
- tamarind_upload_file: Upload files for Tamarind jobs
- tamarind_submit_job: Submit jobs to Tamarind Bio
- task_complete: Call when done with a summary

GENERAL WORKFLOW:
1. Explore the task directory to understand available data
2. Read the task description and requirements carefully
3. Implement the analysis using Python code and/or Tamarind tools
4. For Tamarind tools:
   a. Call tamarind_get_tool_spec to understand required parameters
   b. Upload any required files with tamarind_upload_file
   c. Submit jobs with tamarind_submit_job
5. Save intermediate and final results as JSON/CSV files
6. Call task_complete with a summary when finished

IMPORTANT - TAMARIND TOOL USAGE:
- ALWAYS call tamarind_get_tool_spec before using any tool
- For dict-type parameters (bias_AA_per_residue, omit_AA_per_residue, etc.):
  Pass as actual JSON objects, NOT as stringified JSON.
  CORRECT: "omit_AA_per_residue" should be a dict like {{"A1": "H", "A2": "K"}}
  WRONG:   Do NOT pass as string like "{{...}}" or JSON.stringify
- If a tool call fails with 400 error, check parameter types carefully
- If a tool call fails with 500 error, it may be a transient server issue - the system will retry automatically
- Common tools: esmfold (sequence to structure), proteinmpnn (structure to sequence)

FILE UPLOAD REQUIREMENTS:
- Before using a file (pdbFile, templateFile, etc.) in a Tamarind job, you MUST upload it first
- Call tamarind_upload_file with the file path to upload it
- After upload, use ONLY the filename (not the full path) in job parameters
- Example workflow:
  1. tamarind_upload_file("data/scaffold.pdb")  -> Returns filename "scaffold.pdb"
  2. tamarind_submit_job("proteinmpnn", {{"pdbFile": "scaffold.pdb", ...}})
- Files are session-scoped: only files you uploaded in this session can be used

PYTHON ENVIRONMENT:
- BioPython: PDBParser, ShrakeRupley (for SASA), seq1, SeqIO, NeighborSearch
- NumPy, Pandas, json, Path, os are available
- Use print() to output results
- Current working directory is set to the task directory

Be methodical, save intermediate results, and try alternatives if something fails."""

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": "Please complete the task. Start by exploring the available files."}
    ]
    
    iteration = 0
    task_completed = False
    
    while iteration < max_iterations and not task_completed:
        iteration += 1
        print(f"\n--- Iteration {iteration}/{max_iterations} ---")
        
        try:
            response = client.chat.completions.create(
                model=model,
                messages=messages,
                tools=TOOLS,
                tool_choice="auto"
            )
            
            assistant_message = response.choices[0].message
            messages.append(assistant_message)
            
            if assistant_message.tool_calls:
                for tool_call in assistant_message.tool_calls:
                    tool_name = tool_call.function.name
                    try:
                        arguments = json.loads(tool_call.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {}
                    
                    print(f"  Tool: {tool_name}")
                    if tool_name != "run_python":
                        print(f"  Args: {json.dumps(arguments, indent=2)[:200]}")
                    
                    result = tools.execute_tool(tool_name, arguments)
                    
                    if result.startswith("TASK_COMPLETE:"):
                        task_completed = True
                        print(f"\n{result}")
                    
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": result[:5000]
                    })
                    
                    preview = result[:300] + "..." if len(result) > 300 else result
                    print(f"  Result: {preview}")
            
            elif assistant_message.content:
                print(f"  Assistant: {assistant_message.content[:500]}")
                
        except Exception as e:
            print(f"  Error: {e}")
            traceback.print_exc()
            break
    
    # Save conversation log
    log_path = output_dir / "agent_log.json"
    with open(log_path, "w") as f:
        log_data = [msg.model_dump() if hasattr(msg, "model_dump") else msg for msg in messages]
        json.dump(log_data, f, indent=2, default=str)
    
    # Save session metadata for provenance tracking
    session_info = tools.get_execution_summary()
    session_info["iterations"] = iteration
    session_info["task_completed"] = task_completed
    session_info["model"] = model
    
    session_path = output_dir / "session_info.json"
    with open(session_path, "w") as f:
        json.dump(session_info, f, indent=2, default=str)
    
    print(f"\n{'=' * 60}")
    print(f"Agent finished after {iteration} iterations")
    print(f"Session ID: {session_id}")
    tamarind_info = session_info.get('tamarind', {})
    print(f"Files uploaded: {tamarind_info.get('file_count', 0)}")
    print(f"Jobs submitted: {tamarind_info.get('job_count', 0)}")
    if tamarind_info.get('files_uploaded'):
        print(f"  Files: {tamarind_info['files_uploaded']}")
    if tamarind_info.get('jobs_submitted'):
        print(f"  Jobs: {tamarind_info['jobs_submitted']}")
    print(f"Outputs saved to: {output_dir}")
    print(f"Conversation log: {log_path}")
    print(f"Session info: {session_path}")


def main():
    parser = argparse.ArgumentParser(description="Task-Agnostic AI Agent for Computational Biology")
    parser.add_argument("--task", required=True, help="Task directory name")
    parser.add_argument("--output", help="Output directory name (default: timestamp)")
    parser.add_argument("--max-iterations", type=int, default=100, help="Maximum iterations (default: 20)")
    parser.add_argument("--model", default="gpt-4o", help="OpenAI model (default: gpt-4o)")
    
    args = parser.parse_args()
    run_agent(args.task, args.max_iterations, args.model, args.output)


if __name__ == "__main__":
    main()
