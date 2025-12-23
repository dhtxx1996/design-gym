#!/usr/bin/env python3
"""
Evaluate agent outputs against a rubric.

Usage:
    python eval.py --task ph_sensitive_design --k 3
    python eval.py --task ph_sensitive_design --outputs run1 run2 --k 3
    python eval.py --task ph_sensitive_design --k 5 --parallel  # Parallel across outputs
"""

import os
import re
import json
import asyncio
import argparse
from pathlib import Path
import statistics
from dataclasses import dataclass, field
from typing import Optional

from openai import OpenAI, AsyncOpenAI

# Load .env file and set any missing environment variables
from tamarind_client import _load_env_file
from prompts import build_evaluator_prompt, EVALUATOR_TEMPLATE

_env_vars = _load_env_file()
for _key, _value in _env_vars.items():
    if not os.getenv(_key):
        os.environ[_key] = _value


# ============================================================================
# Tool Failure Analysis
# ============================================================================


@dataclass
class ToolCall:
    """Represents a single tool call."""
    call_id: str
    tool_name: str
    arguments: dict
    result: Optional[str] = None
    is_error: bool = False
    error_type: Optional[str] = None  # "service", "usage", "unknown"


@dataclass 
class ToolFailureAnalysis:
    """Analysis of tool failures in an agent run."""
    total_tool_calls: int = 0
    successful_calls: int = 0
    failed_calls: int = 0  # Any tool call that returned an error/traceback
    
    # Categorized failures (telemetry only; the LLM acts as execution judge)
    service_failures: list = field(default_factory=list)  # Potentially uncontrollable
    usage_failures: list = field(default_factory=list)    # Agent misunderstanding
    recovered_failures: list = field(default_factory=list)  # Failed then succeeded
    
    def to_dict(self) -> dict:
        return {
            "total_tool_calls": self.total_tool_calls,
            "successful_calls": self.successful_calls,
            "failed_calls": self.failed_calls,
            "service_failures": len(self.service_failures),
            "usage_failures": len(self.usage_failures),
            "recovered_failures": len(self.recovered_failures),
        }


def analyze_tool_failures(agent_log_path: Path) -> ToolFailureAnalysis:
    """Analyze agent log for tool calling failures."""
    analysis = ToolFailureAnalysis()
    
    if not agent_log_path.exists():
        return analysis
    
    try:
        log = json.loads(agent_log_path.read_text())
    except:
        return analysis
    
    # Extract tool calls and their results
    tool_calls = {}  # call_id -> ToolCall
    tool_results = {}  # call_id -> result content
    
    for entry in log:
        role = entry.get("role", "")
        
        # Record tool calls from assistant
        if role == "assistant" and entry.get("tool_calls"):
            for tc in entry["tool_calls"]:
                call_id = tc.get("id", "")
                func = tc.get("function", {})
                try:
                    args = json.loads(func.get("arguments", "{}"))
                except:
                    args = {"raw": func.get("arguments", "")}
                
                tool_calls[call_id] = ToolCall(
                    call_id=call_id,
                    tool_name=func.get("name", "unknown"),
                    arguments=args
                )
        
        # Record tool results
        if role == "tool":
            call_id = entry.get("tool_call_id", "")
            content = entry.get("content", "")
            tool_results[call_id] = content
    
    # Match results to calls and analyze
    tool_history = {}  # tool_name -> list of (call_id, success)
    
    for call_id, tc in tool_calls.items():
        result = tool_results.get(call_id, "")
        tc.result = result
        analysis.total_tool_calls += 1
        
        # Check if this is an error
        is_error = any([
            "error" in result.lower()[:200],
            "traceback" in result.lower()[:200],
            "exception" in result.lower()[:200],
            "failed" in result.lower()[:100],
        ])
        
        if is_error:
            tc.is_error = True
            analysis.failed_calls += 1
            # Record generic failure telemetry; the execution judge (LLM)
            # will interpret whether this reflects service issues, usage
            # errors, or other problems based on the raw error text.
            analysis.usage_failures.append({
                "tool": tc.tool_name,
                "call_id": call_id,
                "arguments": tc.arguments,
                "error_snippet": result[:500],
            })
        else:
            analysis.successful_calls += 1
        
        # Track history for recovery detection
        if tc.tool_name not in tool_history:
            tool_history[tc.tool_name] = []
        tool_history[tc.tool_name].append((call_id, not is_error))
    
    # Detect recoveries: tool failed (usage error) then succeeded.
    # This is purely informational telemetry; the execution judge (LLM)
    # decides how to treat these in scoring / force majeure.
    for tool_name, history in tool_history.items():
        had_failure = False
        for i, (call_id, success) in enumerate(history):
            if not success:
                had_failure = True
            elif success and had_failure:
                # Found a recovery: at least one earlier failure for this tool
                # followed by a success.
                analysis.recovered_failures.append(
                    {
                        "tool": tool_name,
                        "description": f"Agent initially failed to use {tool_name} but eventually succeeded",
                    }
                )
                had_failure = False  # Reset for next potential recovery

    return analysis


def parse_rubric(rubric_path: Path) -> tuple[str, dict]:
    """Parse rubric.txt to extract categories and max points."""
    if not rubric_path.exists():
        return "", {}
    
    text = rubric_path.read_text()
    categories = {}
    
    # Try multiple patterns to support different rubric formats
    
    # Pattern 1: "Criterion N: Name" followed by "Max Points: X" on separate line
    pattern1 = r'Criterion\s*\d+:\s*([^\n]+)\n.*?Max Points:\s*(\d+)'
    for match in re.finditer(pattern1, text, re.IGNORECASE | re.DOTALL):
        name = match.group(1).strip()
        points = int(match.group(2))
        key = re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')
        categories[key] = {"name": name, "max": points}
    
    # Pattern 2: "## 1. Category Name (25 pts)" or "## Category (20 points)" (fallback)
    if not categories:
        pattern2 = r'##\s*\d*\.?\s*([^(]+)\s*\((\d+)\s*(?:pts?|points?)\)'
        for match in re.finditer(pattern2, text, re.IGNORECASE):
            name = match.group(1).strip()
            points = int(match.group(2))
            key = re.sub(r'[^a-z0-9]+', '_', name.lower()).strip('_')
            categories[key] = {"name": name, "max": points}
    
    return text, categories


def load_outputs(output_dir: Path) -> dict:
    """Load output files from a directory."""
    outputs = {"path": str(output_dir)}
    
    for f in output_dir.rglob("*.json"):
        if f.stem != "agent_log":
            try:
                outputs[f.stem] = json.loads(f.read_text())
            except:
                pass
    
    for f in output_dir.rglob("*.fa"):
        try:
            outputs[f.stem] = f.read_text()[:2000]
        except:
            pass
    
    # List all files
    outputs["files"] = [str(f.relative_to(output_dir)) for f in output_dir.rglob("*") if f.is_file()]
    
    return outputs


def format_tool_log(output_dir: Path) -> str:
    """Format tool execution log from agent_log.json."""
    agent_log_path = output_dir / "agent_log.json"
    if not agent_log_path.exists():
        return "(No agent execution log found)"
    
    try:
        log_data = json.loads(agent_log_path.read_text())
    except:
        return "(Error reading agent log)"
    
    tool_calls = []
    
    # We need to reconstruct the sequence of calls
    # Filter for assistant tool calls and tool outputs
    
    current_idx = 0
    formatted_log = []
    
    for entry in log_data:
        role = entry.get("role", "")
        
        if role == "assistant" and entry.get("tool_calls"):
            for tc in entry["tool_calls"]:
                func = tc.get("function", {})
                name = func.get("name", "unknown")
                args = func.get("arguments", "{}")
                # Truncate args for display
                if len(args) > 200:
                    args = args[:200] + "...}"
                
                formatted_log.append(f"Step #{current_idx}: {name}({args})")
                current_idx += 1
                
    if not formatted_log:
        return "(No tools called)"
        
    return "\n".join(formatted_log)


def build_failure_context(analysis: ToolFailureAnalysis) -> str:
    """Build context string for the evaluator about tool failures and execution."""
    if analysis.total_tool_calls == 0:
        return ""
    
    lines = []
    
    # Recovered failures - agent should NOT be penalized
    if analysis.recovered_failures:
        lines.append("RECOVERED ERRORS (do NOT penalize - agent learned and fixed these):")
        for rf in analysis.recovered_failures:
            lines.append(f"  - {rf['tool']}: {rf['description']}")
        lines.append("")
    
    # Service failures - potentially uncontrollable (external)
    if analysis.service_failures:
        lines.append("SERVICE-LIKE FAILURES (may indicate external service issues):")
        for sf in analysis.service_failures:
            lines.append(f"  - {sf['tool']}: Service was unavailable/unreachable")
        lines.append("")
    
    # Usage failures that weren't recovered - may factor into evaluation
    unrecovered_usage = [
        uf for uf in analysis.usage_failures
        if not any(rf['tool'] == uf['tool'] for rf in analysis.recovered_failures)
    ]
    if unrecovered_usage:
        lines.append("PERSISTENT TOOL USAGE ERRORS (consider in evaluation):")
        for uf in unrecovered_usage:
            lines.append(f"  - {uf['tool']}: Agent misunderstood tool usage and did not recover")
        lines.append("")
    
    # Summary stats (telemetry for the execution judge)
    lines.append(f"TOOL CALL SUMMARY: {analysis.successful_calls}/{analysis.total_tool_calls} successful")
    
    return "\n".join(lines)


def _build_eval_prompt(task_dir: Path, output_dir: Path, categories: dict, 
                       failure_analysis: Optional[ToolFailureAnalysis] = None) -> tuple[str, int]:
    """Build the evaluation prompt. Returns (prompt, total_max_points)."""
    rubric_text, _ = parse_rubric(task_dir / "rubric.txt")
    question = (task_dir / "question.md").read_text()[:3000] if (task_dir / "question.md").exists() else ""
    outputs = load_outputs(output_dir)
    tool_log_str = format_tool_log(output_dir)
    
    # Build expected JSON format from categories - now with per-category reasoning and navigation links
    if categories:
        category_format = ", ".join(
            f'"{k}": {{"score": <0-{v["max"]}>, "reasoning": "<why this score>", "related_files": ["<filename1>", ...], "related_steps": [<step_idx1>, ...]}}'
            for k, v in categories.items()
        )
        total_max = sum(v["max"] for v in categories.values())
    else:
        category_format = '"overall": {"score": <0-100>, "reasoning": "<explanation>"}'
        total_max = 100
    
    # Build execution/failure context if available
    failure_context = ""
    if failure_analysis:
        failure_context = build_failure_context(failure_analysis)
        if failure_context:
            failure_context = f"\n\nTOOL FAILURE ANALYSIS:\n{failure_context}\n"
    
    # Build prompt from centralized template
    prompt = build_evaluator_prompt(
        question=question,
        rubric_text=rubric_text,
        output_dir_name=output_dir.name,
        files_json=json.dumps(outputs.get('files', []), indent=2),
        contents_json=json.dumps({k: v for k, v in outputs.items() if k not in ['path', 'files']}, indent=2, default=str)[:5000],
        tool_log=tool_log_str,
        failure_context=failure_context,
        category_format=category_format,
        total_max=total_max,
    )
    
    return prompt, total_max


def evaluate_once(task_dir: Path, output_dir: Path, client: OpenAI, model: str, 
                  categories: dict, failure_analysis: Optional[ToolFailureAnalysis] = None) -> dict:
    """Single evaluation round (synchronous)."""
    prompt, _ = _build_eval_prompt(task_dir, output_dir, categories, failure_analysis)
    
    response = client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        response_format={"type": "json_object"}
    )
    
    try:
        return json.loads(response.choices[0].message.content)
    except:
        return {"total": 0, "error": "Parse failed"}


async def evaluate_once_async(task_dir: Path, output_dir: Path, client: AsyncOpenAI, model: str, 
                              categories: dict, failure_analysis: Optional[ToolFailureAnalysis] = None,
                              round_num: int = 0) -> dict:
    """Single evaluation round (async)."""
    prompt, _ = _build_eval_prompt(task_dir, output_dir, categories, failure_analysis)
    
    try:
        response = await client.chat.completions.create(
            model=model,
            messages=[{"role": "user", "content": prompt}],
            response_format={"type": "json_object"}
        )
        result = json.loads(response.choices[0].message.content)
        result["_round"] = round_num
        return result
    except Exception as e:
        return {"total": 0, "error": f"Failed: {e}", "_round": round_num}


def _aggregate_results(all_results: list[dict], categories: dict, 
                       failure_analysis: ToolFailureAnalysis) -> dict:
    """Aggregate results from multiple evaluation rounds."""
    scores = [r.get("total", 0) for r in all_results]
    
    # Average each category - handles {"score": X, "reasoning": "...", "related_files": [...], "related_steps": [...]} format
    category_scores = {}
    for cat in categories:
        cat_data = [r.get(cat, {}) for r in all_results]
        cat_scores = [
            d.get("score", d) if isinstance(d, dict) else d 
            for d in cat_data
        ]
        cat_reasoning = [
            d.get("reasoning", "") if isinstance(d, dict) else ""
            for d in cat_data
        ]
        
        # Aggregate related files and steps across all rounds
        all_related_files = set()
        all_related_steps = set()
        for d in cat_data:
            if isinstance(d, dict):
                files = d.get("related_files")
                if isinstance(files, list):
                    all_related_files.update(files)
                
                steps = d.get("related_steps")
                if isinstance(steps, list):
                    all_related_steps.update(steps)

        category_scores[cat] = {
            "mean": statistics.mean(cat_scores),
            "scores": cat_scores,
            "reasoning": cat_reasoning,
            "related_files": sorted(list(all_related_files)),
            "related_steps": sorted(list(all_related_steps))
        }
    
    # Aggregate execution judge signals
    exec_ok_votes = sum(
        1 for r in all_results
        if isinstance(r.get("execution_ok"), bool) and r.get("execution_ok")
    )
    force_majeure_votes = sum(
        1 for r in all_results
        if isinstance(r.get("force_majeure"), bool) and r.get("force_majeure")
    )
    num_judges = len(all_results)
    exec_ok = exec_ok_votes / max(num_judges, 1) >= 0.5 if num_judges else True
    force_majeure = force_majeure_votes / max(num_judges, 1) >= 0.5 if num_judges else False
    fm_reasons = [
        (r.get("force_majeure_reason") or "").strip()
        for r in all_results
        if isinstance(r.get("force_majeure_reason"), str) and (r.get("force_majeure_reason") or "").strip()
    ]
    seen = set()
    unique_reasons = []
    for reason in fm_reasons:
        if reason not in seen:
            seen.add(reason)
            unique_reasons.append(reason)
    force_majeure_reason = "; ".join(unique_reasons)[:500] if unique_reasons else None

    failure_dict = failure_analysis.to_dict()
    failure_dict.update({
        "execution_judges": num_judges,
        "execution_ok_votes": exec_ok_votes,
        "execution_ok": exec_ok,
        "force_majeure_judges": num_judges,
        "force_majeure_votes": force_majeure_votes,
        "force_majeure": force_majeure,
        "force_majeure_reason": force_majeure_reason,
    })

    return {
        "scores": scores,
        "mean": statistics.mean(scores),
        "std": statistics.stdev(scores) if len(scores) > 1 else 0,
        "breakdown": category_scores,
        "skipped": False,
        "failure_analysis": failure_dict,
    }


def _build_evaluation_entry(eval_result: dict, categories: dict, 
                            eval_model: str, k: int, eval_id: str = None,
                            evaluator_prompt_template: str = None) -> dict:
    """Build an evaluation entry for the manifest."""
    from datetime import datetime
    
    max_points = sum(v["max"] for v in categories.values()) if categories else 100
    mean_score = eval_result.get("mean", 0)
    
    entry = {
        "eval_id": eval_id or f"eval_{datetime.now().strftime('%Y%m%d_%H%M%S')}",
        "evaluated_at": datetime.now().isoformat(),
        "eval_model": eval_model,
        "eval_rounds": k,
        
        "score": {
            "mean": mean_score,
            "std": eval_result.get("std", 0),
            "max": max_points,
            "normalized": mean_score / max_points if max_points > 0 else 0,
            "all_rounds": eval_result.get("scores", []),
        },
        
        "categories": {
            cat_key: {
                "mean": cat_data.get("mean", 0),
                "max": categories.get(cat_key, {}).get("max", 0),
                "scores": cat_data.get("scores", []),
                "reasoning": cat_data.get("reasoning", []),
                "related_files": cat_data.get("related_files", []),
                "related_steps": cat_data.get("related_steps", []),
            }
            for cat_key, cat_data in eval_result.get("breakdown", {}).items()
        },
        
        "failure_analysis": eval_result.get("failure_analysis", {}),
    }
    
    # Save evaluator prompt template for optimization tracking
    if evaluator_prompt_template:
        entry["prompts"] = {
            "evaluator_template": evaluator_prompt_template,
            "evaluator_template_source": "tasks/prompts.py:EVALUATOR_TEMPLATE",
        }
    
    return entry


def _build_rubric_entry(rubric_text: str, categories: dict, 
                        rubric_file_path: str = None) -> dict:
    """Build a rubric entry for the manifest."""
    from datetime import datetime
    
    # Parse criteria with full/partial/zero descriptions
    criteria = []
    for key, cat in categories.items():
        criteria.append({
            "key": key,
            "name": cat.get("name", key),
            "max": cat.get("max", 0),
        })
    
    return {
        "loaded_at": datetime.now().isoformat(),
        "file_path": rubric_file_path,
        "total_points": sum(v["max"] for v in categories.values()) if categories else 100,
        "criteria": criteria,
        "raw_content": rubric_text[:5000],
    }


def save_per_output_eval(output_dir: Path, task_dir: Path, eval_result: dict,
                         categories: dict, eval_model: str, k: int,
                         eval_id: str = None, replace: bool = False) -> Path:
    """Add evaluation to manifest.json in the agent output directory."""
    from datetime import datetime
    
    manifest_path = output_dir / "manifest.json"
    
    # Load existing manifest or create minimal one
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text())
    else:
        # Fallback for old runs without manifest
        manifest = {
            "run_id": output_dir.name,
            "task": task_dir.name,
            "created_at": datetime.now().isoformat(),
            "evaluations": [],
        }
    
    # Build evaluation entry
    rubric_path = task_dir / "rubric.txt"
    rubric_text, _ = parse_rubric(rubric_path)
    rubric_file_path = f"tasks/{task_dir.name}/rubric.txt"
    
    eval_entry = _build_evaluation_entry(
        eval_result, categories, eval_model, k, eval_id,
        evaluator_prompt_template=EVALUATOR_TEMPLATE
    )
    
    # Always update rubric entry with current content
    rubric_entry = _build_rubric_entry(rubric_text, categories, rubric_file_path)
    manifest["rubric"] = rubric_entry
    
    # Add or replace evaluation
    evaluations = manifest.setdefault("evaluations", [])
    if replace and eval_id:
        evaluations = [e for e in evaluations if e.get("eval_id") != eval_id]
        manifest["evaluations"] = evaluations
    evaluations.append(eval_entry)
    
    # Update question if not present
    if not manifest.get("question", {}).get("content"):
        question_path = task_dir / "question.md"
        if question_path.exists():
            manifest["question"] = {
                "content": question_path.read_text()[:5000],
                "file_path": f"tasks/{task_dir.name}/question.md",
            }
    
    # Save updated manifest
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2, default=str)
    
    return manifest_path


async def evaluate_output_async(task_dir: Path, output_dir: Path, client: AsyncOpenAI, 
                                 model: str, categories: dict, k: int,
                                 save_per_output: bool = True,
                                 eval_id: str = None, replace: bool = False) -> tuple[str, dict]:
    """Evaluate a single output directory with k parallel rounds."""
    if not output_dir.exists() or not output_dir.is_dir():
        return output_dir.name, None
    
    # Analyze tool failures
    agent_log_path = output_dir / "agent_log.json"
    failure_analysis = analyze_tool_failures(agent_log_path)
    
    # Run k evaluation rounds in parallel
    tasks = [
        evaluate_once_async(task_dir, output_dir, client, model, categories, failure_analysis, i+1)
        for i in range(k)
    ]
    all_results = await asyncio.gather(*tasks)
    
    # Print results as they complete
    scores = [r.get("total", 0) for r in all_results]
    print(f"  {output_dir.name}: {scores} (mean: {statistics.mean(scores):.1f})")
    
    result = _aggregate_results(all_results, categories, failure_analysis)
    
    # Save to manifest.json
    if save_per_output:
        save_per_output_eval(output_dir, task_dir, result, categories, model, k, eval_id, replace)
    
    return output_dir.name, result


async def evaluate_async(task: str, output_dirs: list[Path], k: int, eval_model: str, 
                         parallel_outputs: bool = False, save_per_output: bool = True,
                         eval_id: str = None, replace: bool = False) -> tuple[dict, dict]:
    """Evaluate outputs with async parallelization.
    
    Args:
        task: Task name
        output_dirs: List of output directories to evaluate
        k: Number of evaluation rounds per output
        eval_model: Model to use for evaluation
        parallel_outputs: If True, evaluate multiple outputs in parallel.
                         If False, only parallelize the k rounds within each output.
        save_per_output: If True, update manifest.json with evaluation.
        eval_id: Optional evaluation ID (auto-generated if not provided).
        replace: If True and eval_id matches existing, replace it.
    """
    client = AsyncOpenAI()
    task_dir = Path(__file__).parent / task
    _, categories = parse_rubric(task_dir / "rubric.txt")
    
    valid_dirs = [d for d in output_dirs if d.exists() and d.is_dir()]
    
    if parallel_outputs:
        # Parallel across all outputs AND all k rounds
        print(f"Running {len(valid_dirs)} outputs × {k} rounds in parallel...")
        tasks = [
            evaluate_output_async(task_dir, output_dir, client, eval_model, categories, k, 
                                  save_per_output, eval_id, replace)
            for output_dir in valid_dirs
        ]
        results_list = await asyncio.gather(*tasks)
        results = {name: result for name, result in results_list if result is not None}
    else:
        # Sequential across outputs, parallel within each output's k rounds
        results = {}
        for output_dir in valid_dirs:
            print(f"\n{output_dir.name}")
            
            # Analyze tool failures
            agent_log_path = output_dir / "agent_log.json"
            failure_analysis = analyze_tool_failures(agent_log_path)
            
            if failure_analysis.total_tool_calls > 0:
                print(f"  Tool calls: {failure_analysis.successful_calls}/{failure_analysis.total_tool_calls} successful")
            
            # Run k rounds in parallel
            tasks = [
                evaluate_once_async(task_dir, output_dir, client, eval_model, categories, failure_analysis, i+1)
                for i in range(k)
            ]
            print(f"  Running {k} evaluation rounds in parallel...")
            all_results = await asyncio.gather(*tasks)
            
            # Sort by round number and print
            all_results_sorted = sorted(all_results, key=lambda r: r.get("_round", 0))
            for r in all_results_sorted:
                print(f"    Round {r.get('_round', '?')}: {r.get('total', 0)}")
            
            result = _aggregate_results(all_results, categories, failure_analysis)
            results[output_dir.name] = result
            
            # Save to manifest.json
            if save_per_output:
                save_per_output_eval(output_dir, task_dir, result, categories, eval_model, k, 
                                     eval_id, replace)
    
    return results, categories


def evaluate(task: str, output_dirs: list[Path], k: int, eval_model: str,
             save_per_output: bool = True, eval_id: str = None, replace: bool = False) -> dict:
    """Evaluate outputs k times (synchronous version).
    
    Args:
        task: Task name
        output_dirs: List of output directories to evaluate
        k: Number of evaluation rounds per output
        eval_model: Model to use for evaluation
        save_per_output: If True, update manifest.json with evaluation.
        eval_id: Optional evaluation ID (auto-generated if not provided).
        replace: If True and eval_id matches existing, replace it.
    """
    client = OpenAI()
    task_dir = Path(__file__).parent / task
    _, categories = parse_rubric(task_dir / "rubric.txt")
    results = {}
    
    for output_dir in output_dirs:
        if not output_dir.exists():
            print(f"Skip: {output_dir} not found")
            continue
        
        if not output_dir.is_dir():
            print(f"Skip: {output_dir} is not a directory")
            continue
        
        print(f"\n{output_dir.name}")
        
        agent_log_path = output_dir / "agent_log.json"
        failure_analysis = analyze_tool_failures(agent_log_path)
        
        if failure_analysis.total_tool_calls > 0:
            print(f"  Tool calls: {failure_analysis.successful_calls}/{failure_analysis.total_tool_calls} successful")
            if failure_analysis.service_failures:
                print(f"  ⚠️  Service failures: {len(failure_analysis.service_failures)}")
            if failure_analysis.recovered_failures:
                print(f"  ✓  Recovered errors: {len(failure_analysis.recovered_failures)}")
            if failure_analysis.usage_failures:
                unrecovered = len(failure_analysis.usage_failures) - len(failure_analysis.recovered_failures)
                if unrecovered > 0:
                    print(f"  ✗  Unrecovered usage errors: {unrecovered}")
        
        scores = []
        all_results = []
        
        for i in range(k):
            result = evaluate_once(task_dir, output_dir, client, eval_model, categories, failure_analysis)
            all_results.append(result)
            score = result.get("total", 0)
            scores.append(score)
            print(f"  Round {i+1}: {score}")
        
        agg_result = _aggregate_results(all_results, categories, failure_analysis)
        results[output_dir.name] = agg_result
        
        # Save to manifest.json
        if save_per_output:
            save_per_output_eval(output_dir, task_dir, agg_result, categories, eval_model, k,
                                 eval_id, replace)
    
    return results, categories


def get_outputs_dir(task: str) -> Path:
    """Get the outputs directory for a task: outputs/{task}/ at repo root."""
    repo_root = Path(__file__).parent.parent
    return repo_root / "outputs" / task


def find_all_output_dirs(task: str) -> list[Path]:
    """Find all output directories for a task."""
    outputs_dir = get_outputs_dir(task)
    if not outputs_dir.exists():
        return []
    
    # Flat structure: outputs/{task}/{run}/
    dirs = []
    for run_dir in sorted(outputs_dir.iterdir()):
        if run_dir.is_dir() and not run_dir.name.startswith('.'):
            if (run_dir / "manifest.json").exists():
                dirs.append(run_dir)
    return dirs


def load_all_manifests(task: str, output_names: list[str] = None) -> list[dict]:
    """Load all manifest.json files from agent outputs for meta-analysis.
    
    Args:
        task: Task name
        output_names: Specific output names to load, or None for all
        
    Returns:
        List of manifest dictionaries, sorted by run_id
    """
    if output_names:
        # Search for specific outputs by name
        outputs_dir = get_outputs_dir(task)
        dirs = []
        for name in output_names:
            if Path(name).is_absolute():
                dirs.append(Path(name))
            else:
                candidate = outputs_dir / name
                if candidate.exists():
                    dirs.append(candidate)
    else:
        dirs = find_all_output_dirs(task)
    
    manifests = []
    for output_dir in dirs:
        manifest_path = output_dir / "manifest.json"
        if manifest_path.exists():
            try:
                manifest = json.loads(manifest_path.read_text())
                manifests.append(manifest)
            except:
                pass
    
    return sorted(manifests, key=lambda m: m.get("run_id", ""))


# Backwards compatibility alias
load_all_eval_results = load_all_manifests


def load_eval_results_as_dataframe(task: str, output_names: list[str] = None, 
                                   expand_evals: bool = True):
    """Load manifests as a pandas DataFrame for analysis.
    
    Args:
        task: Task name
        output_names: Specific output names to load, or None for all
        expand_evals: If True, create one row per evaluation (run × eval).
                      If False, use only the latest evaluation per run.
    
    Returns a flattened DataFrame with:
    - Execution metrics (time, tokens, iterations, etc.)
    - Evaluation scores (mean, std, per-category)
    - Failure analysis metrics
    
    Requires pandas to be installed.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for DataFrame output. Install with: pip install pandas")
    
    manifests = load_all_manifests(task, output_names)
    if not manifests:
        return pd.DataFrame()
    
    rows = []
    for m in manifests:
        execution = m.get("execution", {})
        tokens = execution.get("tokens", {})
        tool_calls = execution.get("tool_calls", {})
        
        # Base row with execution data
        base_row = {
            "task": m.get("task"),
            "run_id": m.get("run_id"),
            "agent_model": execution.get("agent_model"),
            "iterations": execution.get("iterations"),
            "task_completed": execution.get("task_completed"),
            "duration_seconds": execution.get("duration_seconds"),
            "prompt_tokens": tokens.get("prompt"),
            "completion_tokens": tokens.get("completion"),
            "total_tokens": tokens.get("total"),
            "tool_calls_total": tool_calls.get("total"),
            "files_created_count": len(execution.get("files_created", [])),
        }
        
        evaluations = m.get("evaluations", [])
        if not evaluations:
            # No evaluations yet - still include the run
            rows.append(base_row)
            continue
        
        if not expand_evals:
            evaluations = [evaluations[-1]]  # Latest only
        
        for ev in evaluations:
            row = base_row.copy()
            score = ev.get("score", {})
            failure = ev.get("failure_analysis", {})
            
            row.update({
                "eval_id": ev.get("eval_id"),
                "eval_model": ev.get("eval_model"),
                "eval_rounds": ev.get("eval_rounds"),
                "evaluated_at": ev.get("evaluated_at"),
                "score_mean": score.get("mean"),
                "score_std": score.get("std"),
                "score_max": score.get("max"),
                "score_normalized": score.get("normalized"),
                "tool_calls_failed": failure.get("tool_calls_failed", failure.get("failed_calls")),
                "recovered_failures": failure.get("recovered_failures"),
                "execution_ok": failure.get("execution_ok"),
                "force_majeure": failure.get("force_majeure"),
            })
            
            # Add per-category scores
            for cat_key, cat_data in ev.get("categories", {}).items():
                row[f"score_{cat_key}"] = cat_data.get("mean")
            
            rows.append(row)
    
    return pd.DataFrame(rows)


def plot_results(results: dict, output_path: str):
    """Plot scores with error bars."""
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed")
        return
    
    # Filter out skipped results
    valid_results = {k: v for k, v in results.items() if not v.get("skipped", False) and v.get("mean") is not None}
    
    if not valid_results:
        print("No valid results to plot")
        return
    
    names = list(valid_results.keys())
    means = [r["mean"] for r in valid_results.values()]
    stds = [r["std"] for r in valid_results.values()]
    
    fig, ax = plt.subplots(figsize=(max(8, len(names)*0.8), 6))
    bars = ax.bar(range(len(names)), means, yerr=stds, capsize=5, 
                  color='steelblue', edgecolor='black', alpha=0.8)
    
    ax.set_xticks(range(len(names)))
    ax.set_xticklabels(names, rotation=45, ha='right', fontsize=8)
    ax.set_ylabel('Score')
    ax.set_title('Agent Evaluation')
    ax.set_ylim(0, max(means) * 1.2 if means else 100)
    
    for bar, m, s in zip(bars, means, stds):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + s + 2,
                f'{m:.0f}±{s:.0f}', ha='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    print(f"Plot: {output_path}")


def main():
    parser = argparse.ArgumentParser(description="Evaluate agent outputs")
    parser.add_argument("--task", required=True, help="Task name")
    parser.add_argument("--outputs", nargs="+", help="Specific outputs to evaluate")
    parser.add_argument("--k", type=int, default=3, help="Evaluation rounds")
    parser.add_argument("--model", default="gpt-4o", help="Evaluator model")
    parser.add_argument("--output", default="eval_results", help="Output prefix")
    parser.add_argument("--output-dir", type=str, default=None,
                        help="Custom directory to save output files (created if not exists). Default: task directory")
    parser.add_argument("--parallel", action="store_true", 
                        help="Parallelize across outputs (default: only parallelize k rounds per output)")
    parser.add_argument("--sync", action="store_true",
                        help="Use synchronous evaluation (no parallelization, for debugging)")
    parser.add_argument("--no-per-output", action="store_true",
                        help="Don't update manifest.json in each output directory")
    parser.add_argument("--eval-id", type=str, default=None,
                        help="Custom evaluation ID (default: auto-generated timestamp)")
    parser.add_argument("--replace", action="store_true",
                        help="Replace existing evaluation with same eval-id")
    
    args = parser.parse_args()
    task_dir = Path(__file__).parent / args.task
    outputs_dir = get_outputs_dir(args.task)
    save_per_output = not args.no_per_output
    
    # Get output directories
    if args.outputs:
        # Handle specific output names - search in outputs/{task}/{YYYY-MM}/
        output_dirs = []
        for o in args.outputs:
            if Path(o).is_absolute():
                output_dirs.append(Path(o))
            else:
                # Search all month directories for this output name
                found = False
                if outputs_dir.exists():
                    for month_dir in sorted(outputs_dir.iterdir(), reverse=True):
                        if month_dir.is_dir():
                            candidate = month_dir / o
                            if candidate.exists():
                                output_dirs.append(candidate)
                                found = True
                                break
                if not found:
                    print(f"Warning: Output '{o}' not found in {outputs_dir}")
    else:
        output_dirs = find_all_output_dirs(args.task)
    
    if not output_dirs:
        print(f"No outputs found in {outputs_dir}")
        return
    
    print(f"Evaluating {len(output_dirs)} outputs ({args.k} rounds each)")
    if args.eval_id:
        print(f"Evaluation ID: {args.eval_id}")
    if save_per_output:
        print(f"Updating manifest.json in each output directory")
    
    if args.sync:
        print("Using synchronous evaluation (no parallelization)")
        results, categories = evaluate(args.task, output_dirs, args.k, args.model, 
                                        save_per_output, args.eval_id, args.replace)
    else:
        mode = "fully parallel" if args.parallel else "parallel k rounds per output"
        print(f"Using async evaluation ({mode})")
        results, categories = asyncio.run(
            evaluate_async(args.task, output_dirs, args.k, args.model, args.parallel, 
                           save_per_output, args.eval_id, args.replace)
        )
    
    # Determine output directory
    if args.output_dir:
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Saving outputs to: {output_dir}")
    else:
        output_dir = task_dir
    
    # Save
    json_path = output_dir / f"{args.output}.json"
    with open(json_path, "w") as f:
        json.dump({"categories": {k: v["max"] for k, v in categories.items()}, "results": results}, f, indent=2)
    print(f"Saved results: {json_path}")
    
    plot_path = output_dir / f"{args.output}.png"
    plot_results(results, str(plot_path))
    
    # Summary
    cat_names = list(categories.keys())
    header = f"{'Output':<25} {'Total':>10}"
    for c in cat_names:
        header += f" {c[:8]:>8}"
    
    print(f"\n{'='*len(header)}")
    print("RESULTS")
    print("="*len(header))
    print(header)
    print("-"*len(header))
    
    # Sort: non-skipped first by score, then skipped
    def sort_key(item):
        name, r = item
        if r.get("skipped", False) or r.get("mean") is None:
            return (1, 0, name)  # Skipped goes last
        return (0, -r["mean"], name)
    
    for name, r in sorted(results.items(), key=sort_key):
        if r.get("skipped", False) or r.get("mean") is None:
            print(f"{name:<25} {'SKIPPED':>10} (Service unavailable)")
            continue
        
        row = f"{name:<25} {r['mean']:5.1f}±{r['std']:<4.0f}"
        for c in cat_names:
            val = r.get("breakdown", {}).get(c, {}).get("mean", 0)
            row += f" {val:8.1f}"
        print(row)
    
    # Print failure analysis summary
    print(f"\n{'='*60}")
    print("TOOL FAILURE ANALYSIS")
    print("="*60)
    for name, r in results.items():
        fa = r.get("failure_analysis", {})
        if fa.get("total_tool_calls", 0) > 0:
            status = "⛔ SKIPPED" if r.get("skipped") else "✓ EVALUATED"
            print(f"\n{name} [{status}]")
            print(f"  Total calls: {fa.get('total_tool_calls', 0)}, "
                  f"Success: {fa.get('successful_calls', 0)}, "
                  f"Failed: {fa.get('failed_calls', 0)}")
            if fa.get("service_failures", 0) > 0:
                print(f"  ⚠️  Service failures: {fa['service_failures']} (not penalized)")
            if fa.get("recovered_failures", 0) > 0:
                print(f"  ✓  Recovered errors: {fa['recovered_failures']} (not penalized)")
            unrecovered = fa.get("usage_failures", 0) - fa.get("recovered_failures", 0)
            if unrecovered > 0:
                print(f"  ✗  Persistent errors: {unrecovered} (factored into eval)")
    
    # Print per-output file locations and meta-analysis help
    if save_per_output:
        print(f"\n{'='*60}")
        print("MANIFEST FILES (for meta-analysis)")
        print("="*60)
        
        all_tasks_data = {}
        repo_root = Path(__file__).parent.parent
        outputs_root = repo_root / "outputs"
        
        # Update a global index for the dashboard
        for t_dir in outputs_root.iterdir():
            if t_dir.is_dir() and not t_dir.name.startswith('.'):
                task_name = t_dir.name
                run_dirs = find_all_output_dirs(task_name)
                all_tasks_data[task_name] = [rd.name for rd in run_dirs]
        
        index_path = outputs_root / "index.json"
        with open(index_path, "w") as f:
            json.dump(all_tasks_data, f, indent=2)
        print(f"Updated dashboard index: {index_path}")

        for output_dir in output_dirs:
            manifest_path = output_dir / "manifest.json"
            if manifest_path.exists():
                print(f"  {manifest_path}")
        
        print(f"\nTo load all results for analysis:")
        print(f"  from eval import load_all_manifests, load_eval_results_as_dataframe")
        print(f"  manifests = load_all_manifests('{args.task}')")
        print(f"  df = load_eval_results_as_dataframe('{args.task}')")
        print(f"  df = load_eval_results_as_dataframe('{args.task}', expand_evals=True)  # one row per eval")


if __name__ == "__main__":
    main()
