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

from dotenv import load_dotenv
from openai import OpenAI

# Load environment variables
load_dotenv()
load_dotenv(Path(__file__).parent / ".env")

# Add current directory for imports
sys.path.insert(0, str(Path(__file__).parent))
from tamarind_client import TamarindClient
from dag_tracker import DAGTracker


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
    
    def save_dag(self, path: Path = None, show_all_nodes: bool = False) -> None:
        """Save the DAG to JSON and generate visualization.
        
        Args:
            path: Output directory path (default: output_dir)
            show_all_nodes: If True, show all nodes. If False, only show clusters 
                           with at least 3 nodes (default: False)
        """
        if path is None:
            path = self.output_dir
        
        # Save JSON
        json_path = path / "execution_dag.json"
        self.dag.save(json_path)
        
        # Save visualization
        if self.dag.nodes:
            try:
                png_path = path / "execution_dag.png"
                self.dag.visualize(output_path=png_path, title="", show_all=show_all_nodes)
            except Exception as e:
                print(f"Warning: Could not generate DAG visualization: {e}")
    
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
            "tamarind_poll_results": lambda: self.tamarind_poll_results(arguments["job_name"]),
            "task_complete": lambda: f"TASK_COMPLETE: {arguments['summary']}",
        }
        
        # Execute the tool
        result = handlers.get(tool_name, lambda: f"Unknown tool: {tool_name}")()
        
        # Auto-track in DAG (no LLM call needed)
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
    task_dir = tasks_dir / task_name
    name = output_name or datetime.now().strftime("%Y%m%d_%H%M%S")
    output_dir = task_dir / "agent_output" / name
    
    # Generate unique session ID for tracking jobs
    session_id = f"{name}_{datetime.now().strftime('%H%M%S_%f')}"
    
    if not task_dir.exists():
        print(f"Error: Task directory not found: {task_dir}")
        return
    
    # Check for existing output directory
    if output_dir.exists():
        # Check if directory is not empty
        has_contents = any(output_dir.iterdir())
        
        if has_contents:
            if overwrite:
                print(f"Overwriting existing output directory: {output_dir}")
                shutil.rmtree(output_dir)
            else:
                print(f"Error: Output directory already exists and is not empty: {output_dir}")
                print("  Results from previous runs would contaminate evaluation.")
                print("  Options:")
                print("    1. Use --overwrite to delete and replace the existing output")
                print("    2. Use a different --output name")
                print("    3. Manually delete the directory")
                return
    
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
    
    # Load task description from question file
    question_path = Path(question_file) if question_file else task_dir / "question.md"
    if not question_path.is_absolute():
        question_path = task_dir / question_path
    if question_path.exists():
        task_description = question_path.read_text()
        print(f"Question: {question_path.name}")
    else:
        task_description = "Complete the computational biology task in this directory. Explore the files to understand requirements."
        print(f"Question: (default prompt)")
    
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
- tamarind_submit_job: Submit jobs to Tamarind Bio (returns 'pending' if job takes too long)
- tamarind_poll_results: Poll for results of a pending job (use job_name from submit)
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

IF A TAMARIND JOB FAILS:
- Check the error message carefully for hints about what went wrong
- Try simplifying parameters (remove optional ones, use defaults)
- Try a different tool that can achieve the same goal (e.g., esmfold instead of alphafold)
- If bias parameters fail, try without them first to verify the basic job works
- Do NOT give up after one failure - iterate and debug

ASYNC JOB HANDLING (IMPORTANT):
- tamarind_submit_job waits up to 10 minutes for results
- If the job takes longer, it returns status="pending" with the job_name
- When you get a "pending" response:
  1. Note the job_name from the response
  2. Continue with other work (analysis, file prep, etc.)
  3. Later call tamarind_poll_results(job_name) to check if results are ready
  4. If still pending, try again after doing more work
- This allows you to be productive while long-running jobs execute

FILE UPLOAD REQUIREMENTS:
- Before using a file (pdbFile, templateFile, etc.) in a Tamarind job, you MUST upload it first
- First download/create the file, then call tamarind_upload_file to upload it
- After upload, use ONLY the filename (not the full path) in job parameters
- Example workflow:
  1. Download PDB: urllib.request.urlretrieve("https://files.rcsb.org/download/5L33.pdb", "scaffold.pdb")
  2. tamarind_upload_file("scaffold.pdb")  -> Returns filename "scaffold.pdb"
  3. tamarind_submit_job("proteinmpnn", {{"pdbFile": "scaffold.pdb", ...}})
- Files are session-scoped: only files you uploaded in this session can be used

PYTHON ENVIRONMENT:
- Working directory is YOUR session output folder
- OUTPUT_DIR / output_dir: Your session output directory (same as cwd)
- Available: numpy, pandas, BioPython, json, Path, os, urllib

DATA ACCESS:
- Download data from public sources (e.g., RCSB PDB, UniProt) to YOUR session directory
- Working directory is YOUR isolated session folder - all relative paths save there
- Example: urllib.request.urlretrieve("https://files.rcsb.org/download/5L33.pdb", "scaffold.pdb")
- Example: PDBList().retrieve_pdb_file("5L33", file_format="pdb", pdir=output_dir)
- IMPORTANT: All files (downloaded, created, outputs) stay in your session directory

OUTPUT FILES:
- Write to current directory or output_dir (they're the same)
- Do NOT reference any "output/" or "data/" directories
- All your results should be saved via write_file() or in Python to cwd/output_dir

PERSISTENCE - DO NOT GIVE UP EASILY:
- If a computation fails, DEBUG IT - check error messages, fix the code, try again
- If a Tamarind job fails, try different parameters or alternative tools
- If one approach doesn't work, try a different approach to achieve the same goal
- You have many iterations available - use them to iterate and improve
- ONLY call task_complete when you have actually produced the required deliverables
- A failed attempt is NOT a completed task - keep trying until you succeed or exhaust alternatives

Be methodical, save intermediate results, and try alternatives if something fails."""

    messages = [
        {"role": "system", "content": system_prompt},
        {"role": "user", "content": "Please complete the task. Start by exploring the available files."}
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
        print(f"\n--- Iteration {iteration}/{max_iterations} ---")
        
        try:
            response = client.chat.completions.create(
                model=model,
                messages=messages,
                tools=TOOLS,
                tool_choice="auto"
            )
            
            # Track LLM call metrics
            llm_calls += 1
            if response.usage:
                total_prompt_tokens += response.usage.prompt_tokens
                total_completion_tokens += response.usage.completion_tokens
            
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
                    
                    # Track file creation
                    if tool_name == "write_file" and "path" in arguments:
                        files_created.append(arguments["path"])
                    
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
    
    end_time = time.time()
    execution_time_seconds = end_time - start_time
    
    # Save conversation log
    log_path = output_dir / "agent_log.json"
    with open(log_path, "w") as f:
        log_data = [msg.model_dump() if hasattr(msg, "model_dump") else msg for msg in messages]
        json.dump(log_data, f, indent=2, default=str)
    
    # List all files actually created in output directory (excluding logs)
    actual_files = [
        str(f.relative_to(output_dir)) 
        for f in output_dir.rglob("*") 
        if f.is_file() and f.name not in ["agent_log.json", "session_info.json"]
    ]
    
    # Save session metadata for provenance tracking
    session_info = tools.get_execution_summary()
    session_info["iterations"] = iteration
    session_info["task_completed"] = task_completed
    session_info["model"] = model
    
    # Add execution metrics
    session_info["execution_metrics"] = {
        "execution_time_seconds": round(execution_time_seconds, 2),
        "execution_time_formatted": f"{int(execution_time_seconds // 60)}m {int(execution_time_seconds % 60)}s",
        "llm_calls": llm_calls,
        "total_prompt_tokens": total_prompt_tokens,
        "total_completion_tokens": total_completion_tokens,
        "total_tokens": total_prompt_tokens + total_completion_tokens,
        "files_created": actual_files,
        "files_created_count": len(actual_files),
    }
    
    # Add DAG metadata
    session_info["execution_dag"] = tools.dag.to_dict().get("metadata", {})
    
    session_path = output_dir / "session_info.json"
    with open(session_path, "w") as f:
        json.dump(session_info, f, indent=2, default=str)
    
    # Save execution DAG
    tools.save_dag(show_all_nodes=show_all_dag_nodes)
    dag_stats = tools.dag.to_dict().get("metadata", {})
    
    print(f"\n{'=' * 60}")
    print(f"Agent finished after {iteration} iterations")
    print(f"Execution time: {session_info['execution_metrics']['execution_time_formatted']}")
    print(f"LLM calls: {llm_calls}")
    print(f"Tokens: {total_prompt_tokens:,} prompt + {total_completion_tokens:,} completion = {total_prompt_tokens + total_completion_tokens:,} total")
    print(f"Files created: {len(actual_files)}")
    print(f"Session ID: {session_id}")
    tamarind_info = session_info.get('tamarind', {})
    print(f"Tamarind files uploaded: {tamarind_info.get('file_count', 0)}")
    print(f"Tamarind jobs submitted: {tamarind_info.get('job_count', 0)}")
    if tamarind_info.get('files_uploaded'):
        print(f"  Files: {tamarind_info['files_uploaded']}")
    if tamarind_info.get('jobs_submitted'):
        print(f"  Jobs: {tamarind_info['jobs_submitted']}")
    print(f"Execution DAG: {dag_stats.get('node_count', 0)} nodes, {dag_stats.get('edge_count', 0)} edges")
    print(f"Outputs saved to: {output_dir}")
    print(f"Conversation log: {log_path}")
    print(f"Session info: {session_path}")
    if dag_stats.get('node_count', 0) > 0:
        print(f"DAG visualization: {output_dir / 'execution_dag.png'}")


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
