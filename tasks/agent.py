#!/usr/bin/env python3
"""
Task-Agnostic AI Agent for Computational Biology

This agent uses OpenAI's function calling to work on computational biology
tasks defined in task directories. Each task should have:
- question.md: Task description and requirements

Usage:
    python agent.py --task <task_name> [--max-iterations 20]
"""

import os
import sys
import json
import shutil
import argparse
import traceback
from pathlib import Path
from datetime import datetime
from typing import Optional

# Add current directory for imports
sys.path.insert(0, str(Path(__file__).parent))

# Load environment variables from .env file (fallback if not in environment)
from tamarind_client import _load_env_file, TamarindClient
from prompts import build_agent_system_prompt, AGENT_INITIAL_MESSAGE

# Load .env and set any missing environment variables
_env_vars = _load_env_file()
for _key, _value in _env_vars.items():
    if not os.getenv(_key):
        os.environ[_key] = _value

from openai import OpenAI
from dag_tracker import DAGTracker


# =============================================================================
# Color Utilities for Terminal Output
# =============================================================================

class Colors:
    """ANSI color codes for terminal output."""
    # Reset
    RESET = '\033[0m'
    
    # Regular colors
    BLACK = '\033[30m'
    RED = '\033[31m'
    GREEN = '\033[32m'
    YELLOW = '\033[33m'
    BLUE = '\033[34m'
    MAGENTA = '\033[35m'
    CYAN = '\033[36m'
    WHITE = '\033[37m'
    
    # Bright colors
    BRIGHT_BLACK = '\033[90m'
    BRIGHT_RED = '\033[91m'
    BRIGHT_GREEN = '\033[92m'
    BRIGHT_YELLOW = '\033[93m'
    BRIGHT_BLUE = '\033[94m'
    BRIGHT_MAGENTA = '\033[95m'
    BRIGHT_CYAN = '\033[96m'
    BRIGHT_WHITE = '\033[97m'
    
    # Styles
    BOLD = '\033[1m'
    DIM = '\033[2m'
    UNDERLINE = '\033[4m'
    
    @staticmethod
    def colorize(text: str, color: str, bold: bool = False) -> str:
        """Apply color and optional bold styling to text."""
        style = Colors.BOLD if bold else ""
        return f"{style}{color}{text}{Colors.RESET}"
    
    @staticmethod
    def header(text: str) -> str:
        """Format header text."""
        return Colors.colorize(text, Colors.BRIGHT_CYAN, bold=True)
    
    @staticmethod
    def success(text: str) -> str:
        """Format success message."""
        return Colors.colorize(text, Colors.BRIGHT_GREEN)
    
    @staticmethod
    def error(text: str) -> str:
        """Format error message."""
        return Colors.colorize(text, Colors.BRIGHT_RED, bold=True)
    
    @staticmethod
    def warning(text: str) -> str:
        """Format warning message."""
        return Colors.colorize(text, Colors.BRIGHT_YELLOW)
    
    @staticmethod
    def info(text: str) -> str:
        """Format info message."""
        return Colors.colorize(text, Colors.BRIGHT_BLUE)
    
    @staticmethod
    def tool(text: str) -> str:
        """Format tool call."""
        return Colors.colorize(text, Colors.BRIGHT_YELLOW)
    
    @staticmethod
    def iteration(text: str) -> str:
        """Format iteration marker."""
        return Colors.colorize(text, Colors.BRIGHT_MAGENTA, bold=True)
    
    @staticmethod
    def highlight(text: str) -> str:
        """Highlight important text."""
        return Colors.colorize(text, Colors.BRIGHT_CYAN)
    
    @staticmethod
    def dim(text: str) -> str:
        """Dim less important text."""
        return Colors.colorize(text, Colors.DIM)


# =============================================================================
# Tool Definitions for OpenAI Function Calling
# =============================================================================

TOOLS = [
    {
        "type": "function",
        "function": {
            "name": "read_file",
            "description": "Read the contents of a file from your session output directory.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Filename or relative path within your session output directory."
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
            "description": "Write content to a file in YOUR session output directory. All outputs are isolated per run.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Filename or relative path for output (e.g., 'results.json', 'analysis/data.csv')"
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
            "description": "List files in a directory. Use 'data' for input files, '' for task root.",
            "parameters": {
                "type": "object",
                "properties": {
                    "path": {
                        "type": "string",
                        "description": "Directory path: 'data' for inputs, '' for task root, or path in your outputs"
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
            "description": "Execute Python code. CRITICAL: (1) Each call is ISOLATED - variables do NOT persist between calls. (2) You MUST use print() to see output. (3) Save results to JSON files if needed in later calls. Available: numpy, pandas, BioPython (PDBParser, ShrakeRupley, seq1, SeqIO, PDBList), json, Path, os, urllib.",
            "parameters": {
                "type": "object",
                "properties": {
                    "code": {
                        "type": "string",
                        "description": "Python code to execute. Use print() for output. Variables don't persist - save to files."
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
            "description": "Submit a job to Tamarind Bio and wait for results. Pass all parameters inside the 'params' dict. If the job takes too long, it returns a 'pending' status with the job_name - use tamarind_poll_results to check later.",
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
            "name": "tamarind_poll_results",
            "description": "Poll for results of a previously submitted Tamarind job. Use this when tamarind_submit_job returned a 'pending' status. Can be called multiple times until job completes.",
            "parameters": {
                "type": "object",
                "properties": {
                    "job_name": {
                        "type": "string",
                        "description": "The job_name returned from tamarind_submit_job"
                    }
                },
                "required": ["job_name"]
            }
        }
    },
    {
        "type": "function",
        "function": {
            "name": "task_complete",
            "description": "Call ONLY when you have successfully completed the task AND produced the required deliverables. Do NOT call if you encountered failures you haven't resolved - keep trying alternatives instead.",
            "parameters": {
                "type": "object",
                "properties": {
                    "summary": {
                        "type": "string",
                        "description": "Summary of completed work, deliverables produced, and key findings. List what was achieved, not what failed."
                    }
                },
                "required": ["summary"]
            }
        }
    },
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
        
        # DAG tracking for execution provenance
        self.dag = DAGTracker()
        
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
        # Only read from agent's own output directory
        full_path = self.output_dir / path
        if not full_path.exists():
            return f"Error: File not found: {path}. Download data first using urllib or PDBList in run_python."
        
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
        # Only list within agent's own output directory
        full_path = self.output_dir / path if path else self.output_dir
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
        # Set working directory to agent's output dir so relative writes go there
        os.chdir(self.output_dir)
        
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
            from Bio.PDB import PDBParser, NeighborSearch, PDBIO, Selection, PDBList
            from Bio.PDB.SASA import ShrakeRupley
            from Bio.SeqUtils import seq1
            from Bio import SeqIO
            exec_globals.update({
                "Bio": Bio, "PDBParser": PDBParser, "ShrakeRupley": ShrakeRupley,
                "NeighborSearch": NeighborSearch, "PDBIO": PDBIO, "Selection": Selection,
                "PDBList": PDBList, "seq1": seq1, "SeqIO": SeqIO
            })
        except ImportError:
            pass
        
        # Provide output path - working directory is set to OUTPUT_DIR
        exec_globals["OUTPUT_DIR"] = self.output_dir
        exec_globals["output_dir"] = str(self.output_dir)
        
        # Add urllib for downloading data
        import urllib.request
        exec_globals["urllib"] = urllib
        
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
        # Only upload files from agent's output directory
        full_path = self.output_dir / filepath
        if not full_path.exists():
            return f"Error: File not found: {filepath}. Download data first using urllib or PDBList in run_python."
        
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
            
            # Check if job is still pending (didn't complete within timeout)
            final_status = result.get("final_status", {})
            if final_status.get("status") == "pending":
                # Job is still running - return info for async polling
                return json.dumps({
                    "status": "pending",
                    "job_name": result.get("job_name"),
                    "message": (
                        f"Job '{result.get('job_name')}' was submitted but is still running. "
                        f"Use tamarind_poll_results with job_name='{result.get('job_name')}' to check for results later. "
                        f"You can continue with other work in the meantime."
                    ),
                    "tool": tool_name,
                    "submitted_at": result.get("submitted_at"),
                }, indent=2, default=str)
            
            # Job completed - try to download results
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
    
    def tamarind_poll_results(self, job_name: str) -> str:
        """Poll for results of a previously submitted job."""
        try:
            result = self.tamarind.poll_job_results(job_name, timeout=120, poll_interval=10)
            
            if result.get("status") == "completed":
                # Job completed - try to download results
                try:
                    download_path = self.tamarind.download_results(
                        job_name, 
                        output_dir=str(self.output_dir / "tamarind_results")
                    )
                    result["downloaded_to"] = str(download_path)
                    result["message"] = f"Job completed! Results downloaded to {download_path}"
                except Exception as download_error:
                    result["download_error"] = str(download_error)
                    result["note"] = "Job completed but result download failed. Results may be available on Tamarind web UI."
                
                return json.dumps(result, indent=2, default=str)
            
            elif result.get("status") == "pending":
                # Still running
                return json.dumps({
                    "status": "pending",
                    "job_name": job_name,
                    "message": (
                        f"Job '{job_name}' is still running. "
                        f"Try again later with tamarind_poll_results, or continue with other work."
                    ),
                    "elapsed_seconds": result.get("elapsed_seconds"),
                }, indent=2, default=str)
            
            else:
                # Failed or other status
                return json.dumps(result, indent=2, default=str)
                
        except Exception as e:
            return f"Error polling results: {e}"
    
    # =========================================================================
    # Automatic DAG Tracking
    # =========================================================================
    
    def _auto_track_tool(self, tool_name: str, arguments: dict, result: str) -> None:
        """Automatically track tool call in the DAG."""
        try:
            # Generate unique function ID
            func_count = len([n for n in self.dag.nodes.values() if n.type == "function"])
            func_id = f"{tool_name}_{func_count + 1}"
            
            # Skip tracking for some meta tools
            if tool_name in ("list_directory", "task_complete"):
                return
            
            # Helper to get just filename from path
            def basename(p: str) -> str:
                return Path(p).name
            
            # Determine inputs and outputs based on tool type
            inputs = []
            outputs = []
            
            if tool_name == "read_file":
                path = arguments.get("path", "file")
                artifact_id = f"file:{path}"
                if artifact_id not in self.dag.nodes:
                    self.dag.add_artifact(basename(path), artifact_id=artifact_id)
                inputs.append(artifact_id)
                
            elif tool_name == "write_file":
                path = arguments.get("path", "file")
                artifact_id = f"file:{path}"
                if artifact_id not in self.dag.nodes:
                    self.dag.add_artifact(basename(path), artifact_id=artifact_id)
                outputs.append(artifact_id)
                
            elif tool_name == "run_python":
                # Track code execution - look for file operations in code
                code = arguments.get("code", "")
                # Check for file writes in code
                import re
                writes = re.findall(r'\.write_text\(["\']([^"\']+)', code)
                writes += re.findall(r'urlretrieve\([^,]+,\s*["\']([^"\']+)', code)
                writes += re.findall(r'open\(["\']([^"\']+)["\'],\s*["\']w', code)
                writes += re.findall(r'to_json\(["\']([^"\']+)', code)
                writes += re.findall(r'json\.dump\([^,]+,\s*open\(["\']([^"\']+)', code)
                for path in writes:
                    artifact_id = f"file:{path}"
                    if artifact_id not in self.dag.nodes:
                        self.dag.add_artifact(basename(path), artifact_id=artifact_id)
                    outputs.append(artifact_id)
                    
            elif tool_name == "tamarind_upload_file":
                path = arguments.get("filepath", "file")
                artifact_id = f"file:{path}"
                if artifact_id not in self.dag.nodes:
                    self.dag.add_artifact(basename(path), artifact_id=artifact_id)
                inputs.append(artifact_id)
                # Output is the uploaded reference
                upload_id = f"tamarind:{path}"
                self.dag.add_artifact(basename(path), artifact_id=upload_id)
                outputs.append(upload_id)
                
            elif tool_name == "tamarind_submit_job":
                tool = arguments.get("tool_name", "job")
                func_id = f"{tool}_{func_count + 1}"
                # Input is any uploaded file referenced
                params = arguments.get("params", {})
                for key, val in params.items():
                    if isinstance(val, str) and val.endswith(('.pdb', '.fa', '.fasta', '.json')):
                        upload_id = f"tamarind:{val}"
                        if upload_id in self.dag.nodes:
                            inputs.append(upload_id)
                # Output is results
                if "job_name" in result:
                    result_id = f"results:{tool}"
                    self.dag.add_artifact(f"{tool}_results", artifact_id=result_id)
                    outputs.append(result_id)
                    
            elif tool_name == "tamarind_poll_results":
                job_name = arguments.get("job_name", "")
                tool = job_name.split("_")[0] if "_" in job_name else "job"
                result_id = f"results:{tool}"
                if result_id in self.dag.nodes:
                    inputs.append(result_id)
                # Downloaded files as outputs
                if "downloaded_to" in result:
                    dl_id = f"downloaded:{tool}"
                    self.dag.add_artifact(f"{tool}_downloaded", artifact_id=dl_id)
                    outputs.append(dl_id)
            
            # Add function node if we have any i/o
            if inputs or outputs:
                self.dag.add_function(tool_name, function_id=func_id)
                for inp in inputs:
                    try:
                        self.dag.add_edge(inp, func_id)
                    except:
                        pass
                for out in outputs:
                    try:
                        self.dag.add_edge(func_id, out)
                    except:
                        pass
                        
        except Exception as e:
            # Don't fail the tool call if DAG tracking fails
            print(f"[dag] Warning: tracking failed: {e}")
    
    def execute_tool(self, tool_name: str, arguments: dict, reasoning: str = None, tokens: dict = None) -> str:
        import time as _time
        handlers = {
            "read_file": lambda: self.read_file(arguments["path"]),
            "write_file": lambda: self.write_file(arguments["path"], arguments["content"]),
            "list_directory": lambda: self.list_directory(arguments["path"]),
            "run_python": lambda: self.run_python(arguments["code"]),
            "tamarind_list_tools": lambda: self.tamarind_list_tools(),
            "tamarind_get_tool_spec": lambda: self.tamarind_get_tool_spec(arguments["tool_name"]),
            "tamarind_upload_file": lambda: self.tamarind_upload_file(arguments["filepath"]),
            "tamarind_submit_job": lambda: self._handle_tamarind_submit(arguments),
            "tamarind_poll_results": lambda: self.tamarind_poll_results(arguments["job_name"]),
            "task_complete": lambda: f"TASK_COMPLETE: {arguments['summary']}",
        }
        
        # Execute and time the tool
        start = _time.time()
        result = handlers.get(tool_name, lambda: f"Unknown tool: {tool_name}")()
        duration_ms = int((_time.time() - start) * 1000)
        
        # Track tool call for manifest
        is_error = "error" in result.lower()[:200] or "traceback" in result.lower()[:200]
        self.tool_calls.append({
            "idx": len(self.tool_calls),
            "name": tool_name,
            "arguments": arguments,
            "reasoning": reasoning,
            "tokens": tokens,
            "timestamp": datetime.now().isoformat(),
            "duration_ms": duration_ms,
            "status": "error" if is_error else "success",
            "preview": result[:150].replace("\n", " "),
            "result": result if len(result) < 10000 else result[:10000] + "\n... [truncated]"
        })
        
        # Auto-track in DAG
        self._auto_track_tool(tool_name, arguments, result)
        
        return result
    
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

def run_agent(task_name: str, max_iterations: int = 20, model: str = "gpt-4o", output_name: str = None, overwrite: bool = False, question_file: str = None, show_all_dag_nodes: bool = False):
    """Run the agent on a task."""
    tasks_dir = Path(__file__).parent
    repo_root = tasks_dir.parent
    task_dir = tasks_dir / task_name
    
    # Build output path: outputs/{task}/{YYYYMMDD_HHMMSS}_{name}/
    now = datetime.now()
    timestamp = now.strftime("%Y%m%d_%H%M%S")
    run_name = f"{timestamp}_{output_name}" if output_name else timestamp
    output_dir = repo_root / "outputs" / task_name / run_name
    
    # Generate unique session ID for tracking jobs
    session_id = f"{run_name}_{now.strftime('%H%M%S_%f')}"
    
    if not task_dir.exists():
        print(f"{Colors.error('Error:')} Task directory not found: {task_dir}")
        return
    
    # Check for existing output directory
    if output_dir.exists():
        # Check if directory is not empty
        has_contents = any(output_dir.iterdir())
        
        if has_contents:
            if overwrite:
                print(f"{Colors.warning('Overwriting existing output directory:')} {output_dir}")
                shutil.rmtree(output_dir)
            else:
                print(f"{Colors.error('Error:')} Output directory already exists and is not empty: {output_dir}")
                print(f"  {Colors.warning('Results from previous runs would contaminate evaluation.')}")
                print(f"  {Colors.info('Options:')}")
                print(f"    1. Use --overwrite to delete and replace the existing output")
                print(f"    2. Use a different --output name")
                print(f"    3. Manually delete the directory")
                return
    
    print(Colors.header("=" * 60))
    print(Colors.header("AI Agent for Computational Biology"))
    print(Colors.header("=" * 60))
    print(f"{Colors.info('Task:')} {task_name}")
    print(f"{Colors.info('Output:')} {output_dir}")
    print(f"{Colors.info('Session ID:')} {session_id}")
    print(f"{Colors.info('Model:')} {model}")
    print(f"{Colors.info('Max iterations:')} {max_iterations}")
    print()
    
    client = OpenAI()
    tools = AgentTools(task_dir, output_dir, session_id)
    
    # Load task description from question file
    question_path = Path(question_file) if question_file else task_dir / "question.md"
    if not question_path.is_absolute():
        question_path = task_dir / question_path
    if question_path.exists():
        task_description = question_path.read_text()
        print(f"{Colors.info('Question:')} {question_path.name}")
    else:
        task_description = "Complete the computational biology task in this directory. Explore the files to understand requirements."
        print(f"{Colors.info('Question:')} (default prompt)")
    
    # Task-agnostic system prompt
    system_prompt = build_agent_system_prompt(task_description)

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": AGENT_INITIAL_MESSAGE}
    ]
    
    iteration = 0
    task_completed = False
    
    # Track execution metrics
    import time
    start_time = time.time()
    llm_calls = 0
    total_prompt_tokens = 0
    total_completion_tokens = 0
    files_created = []
    
    while iteration < max_iterations and not task_completed:
        iteration += 1
        print(f"\n{Colors.iteration(f'--- Iteration {iteration}/{max_iterations} ---')}")
        
        try:
            response = client.chat.completions.create(
                model=model,
                messages=messages,
                tools=TOOLS,
                tool_choice="auto"
            )
            
            # Track LLM call metrics
            llm_calls += 1
            turn_tokens = None
            if response.usage:
                turn_tokens = {
                    "prompt": response.usage.prompt_tokens,
                    "completion": response.usage.completion_tokens,
                    "total": response.usage.total_tokens
                }
                total_prompt_tokens += response.usage.prompt_tokens
                total_completion_tokens += response.usage.completion_tokens
            
            assistant_message = response.choices[0].message
            messages.append(assistant_message)
            
            # Display assistant's reasoning if present
            reasoning = assistant_message.content
            if reasoning:
                print(f"\n  {Colors.info('Assistant Reasoning:')}")
                for line in reasoning.split('\n'):
                    print(f"    {line}")
            
            if assistant_message.tool_calls:
                for tool_call in assistant_message.tool_calls:
                    tool_name = tool_call.function.name
                    try:
                        arguments = json.loads(tool_call.function.arguments)
                    except json.JSONDecodeError:
                        arguments = {}
                    
                    print(f"  {Colors.tool('Tool:')} {Colors.highlight(tool_name)}")
                    if tool_name == "run_python":
                        code = arguments.get("code", "")
                        print(f"  {Colors.dim('Code:')}")
                        for line in code.split('\n'):
                            print(f"    {Colors.colorize(line, Colors.GREEN)}")
                    else:
                        print(f"  {Colors.dim('Args:')} {Colors.dim(json.dumps(arguments, indent=2)[:200])}")
                    
                    # Track file creation
                    if tool_name == "write_file" and "path" in arguments:
                        files_created.append(arguments["path"])
                    
                    # Pass reasoning and tokens to execute_tool for observability
                    result = tools.execute_tool(
                        tool_name, 
                        arguments, 
                        reasoning=reasoning if tool_call == assistant_message.tool_calls[0] else None,
                        tokens=turn_tokens if tool_call == assistant_message.tool_calls[0] else None
                    )
                    
                    if result.startswith("TASK_COMPLETE:"):
                        task_completed = True
                        print(f"\n{Colors.success(result)}")
                    
                    messages.append({
                        "role": "tool",
                        "tool_call_id": tool_call.id,
                        "content": result[:5000]
                    })
                    
                    preview = result[:300] + "..." if len(result) > 300 else result
                    print(f"  {Colors.dim('Result:')} {Colors.dim(preview)}")
                
        except Exception as e:
            print(f"  {Colors.error('Error:')} {e}")
            traceback.print_exc()
            break
    
    end_time = time.time()
    execution_time_seconds = end_time - start_time
    
    # Save conversation log
    log_path = output_dir / "agent_log.json"
    with open(log_path, "w") as f:
        log_data = [msg.model_dump() if hasattr(msg, "model_dump") else msg for msg in messages]
        json.dump(log_data, f, indent=2, default=str)
    
    # List all files created (excluding meta files)
    meta_files = {"agent_log.json", "manifest.json", "dag.png"}
    actual_files = [
        str(f.relative_to(output_dir)) 
        for f in output_dir.rglob("*") 
        if f.is_file() and f.name not in meta_files
    ]
    
    # Build unified manifest.json
    tamarind_summary = tools.tamarind.get_session_summary() if tools._tamarind else {}
    dag_data = tools.dag.to_dict()
    
    # Collect system environment info
    import platform
    system_info = {
        "os": platform.system(),
        "os_release": platform.release(),
        "python_version": platform.python_version(),
        "libraries": {}
    }
    for lib in ["numpy", "pandas", "Bio", "openai"]:
        try:
            import importlib
            m = importlib.import_module(lib)
            system_info["libraries"][lib] = getattr(m, "__version__", "unknown")
        except ImportError:
            system_info["libraries"][lib] = "not_installed"

    manifest = {
        "run_id": run_name,
        "task": task_name,
        "created_at": datetime.now().isoformat(),
        
        # Execution data
        "execution": {
            "session_id": session_id,
            "agent_model": model,
            "iterations": iteration,
            "task_completed": task_completed,
            "start_time": tools.execution_start.isoformat(),
            "end_time": datetime.now().isoformat(),
            "duration_seconds": round(execution_time_seconds, 2),
            "tokens": {
                "prompt": total_prompt_tokens,
                "completion": total_completion_tokens,
                "total": total_prompt_tokens + total_completion_tokens,
            },
            "tool_calls": {
                "total": len(tools.tool_calls),
                "by_name": {},
            },
            "tamarind": tamarind_summary,
            "files_created": actual_files,
            "system": system_info,
        },
        
        # DAG (embedded)
        "dag": {
            "nodes": dag_data.get("nodes", []),
            "edges": dag_data.get("edges", []),
        },
        
        # Tool call log
        "tool_log": tools.tool_calls,
        
        # Prompts used (for optimization tracking)
        "prompts": {
            "agent_system": system_prompt,
            "agent_system_source": "tasks/agent.py:run_agent",
            "agent_initial": "Please complete the task. Start by exploring the available files.",
        },
        
        # Question/task input
        "question": {
            "content": task_description,
            "file_path": str(question_path.relative_to(repo_root)) if question_path.exists() else None,
        },
        
        # Placeholders for eval.py
        "evaluations": [],
        "rubrics": [],
    }
    
    # Count tool calls by name
    for tc in tools.tool_calls:
        name_key = tc.get("name", "unknown")
        manifest["execution"]["tool_calls"]["by_name"][name_key] = \
            manifest["execution"]["tool_calls"]["by_name"].get(name_key, 0) + 1
    
    manifest_path = output_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    
    # Save DAG visualization only
    dag_stats = dag_data.get("metadata", {})
    if tools.dag.nodes:
        try:
            png_path = output_dir / "dag.png"
            tools.dag.visualize(output_path=png_path, title="", show_all=show_all_dag_nodes)
        except Exception as e:
            print(f"Warning: Could not generate DAG visualization: {e}")
    
    # Print summary
    exec_time_fmt = f"{int(execution_time_seconds // 60)}m {int(execution_time_seconds % 60)}s"
    print(f"\n{Colors.header('=' * 60)}")
    print(Colors.header("Agent finished"))
    print(f"{Colors.highlight('Iterations:')} {iteration}")
    print(f"{Colors.highlight('Execution time:')} {exec_time_fmt}")
    print(f"{Colors.highlight('LLM calls:')} {llm_calls}")
    print(f"{Colors.highlight('Tokens:')} {total_prompt_tokens:,} prompt + {total_completion_tokens:,} completion = {total_prompt_tokens + total_completion_tokens:,} total")
    print(f"{Colors.highlight('Files created:')} {len(actual_files)}")
    print(f"{Colors.info('Session ID:')} {session_id}")
    print(f"{Colors.info('Tamarind files uploaded:')} {tamarind_summary.get('file_count', 0)}")
    print(f"{Colors.info('Tamarind jobs submitted:')} {tamarind_summary.get('job_count', 0)}")
    print(f"{Colors.info('Execution DAG:')} {dag_stats.get('node_count', 0)} nodes, {dag_stats.get('edge_count', 0)} edges")
    print(f"{Colors.info('Outputs saved to:')} {output_dir}")
    print(f"{Colors.dim('Manifest:')} {manifest_path}")
    print(f"{Colors.dim('Conversation log:')} {log_path}")
    if dag_stats.get('node_count', 0) > 0:
        print(f"{Colors.dim('DAG visualization:')} {output_dir / 'dag.png'}")


def main():
    parser = argparse.ArgumentParser(description="Task-Agnostic AI Agent for Computational Biology")
    parser.add_argument("--task", required=True, help="Task directory name")
    parser.add_argument("--question", help="Question file path (default: question.md)")
    parser.add_argument("--output", help="Output directory name (default: timestamp)")
    parser.add_argument("--max-iterations", type=int, default=100, help="Maximum iterations (default: 100)")
    parser.add_argument("--model", default="gpt-4o", help="OpenAI model (default: gpt-4o)")
    parser.add_argument("--overwrite", action="store_true", help="Delete and overwrite existing output directory")
    parser.add_argument("--show-all-dag-nodes", action="store_true", 
                        help="Show all nodes in DAG visualization (default: only clusters with 3+ nodes)")
    
    args = parser.parse_args()
    run_agent(args.task, args.max_iterations, args.model, args.output, args.overwrite, args.question, args.show_all_dag_nodes)


if __name__ == "__main__":
    main()
