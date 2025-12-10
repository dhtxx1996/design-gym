#!/usr/bin/env python3
"""
Evaluate agent outputs against a rubric.

Usage:
    python eval.py --task ph_sensitive_design --k 3
    python eval.py --task ph_sensitive_design --outputs run1 run2 --k 3
    python eval.py --task ph_sensitive_design --k 5 --parallel  # Parallel across outputs
"""

import re
import json
import asyncio
import argparse
from pathlib import Path
import statistics
from dataclasses import dataclass, field
from typing import Optional

from dotenv import load_dotenv
from openai import OpenAI, AsyncOpenAI

load_dotenv()
load_dotenv(Path(__file__).parent / ".env")


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
    
    # Build expected JSON format from categories - now with per-category reasoning
    if categories:
        category_format = ", ".join(
            f'"{k}": {{"score": <0-{v["max"]}>, "reasoning": "<why this score>"}}'
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
    
    prompt = f"""Evaluate this AI agent's solution to a computational biology task.

TASK:
{question}

RUBRIC:
{rubric_text}

AGENT OUTPUT DIRECTORY: {output_dir.name}

FILES PRODUCED:
{json.dumps(outputs.get('files', []), indent=2)}

FILE CONTENTS:
{json.dumps({k: v for k, v in outputs.items() if k not in ['path', 'files']}, indent=2, default=str)[:5000]}
{failure_context}
EVALUATION GUIDELINES FOR TOOL FAILURES AND EXECUTION CONDITIONS:
1. Use the TOOL FAILURE ANALYSIS section as factual telemetry about tool usage and errors.
2. If tools/services were clearly unavailable (e.g., repeated network/HTTP 5xx errors), do NOT penalize
   the agent for those "force majeure" conditions.
3. RECOVERED ERRORS: If the agent initially misunderstood a tool but eventually figured it out and 
   succeeded, do NOT penalize for the initial failures. Learning and recovery is expected behavior.
4. PERSISTENT MISUSE: If the agent consistently failed to use a tool correctly despite multiple 
   attempts (and never recovered), this MAY factor into the evaluation - but focus on whether the 
   agent achieved the task goals through alternative means.
5. FOCUS ON OUTCOMES: Primarily evaluate based on the final outputs and whether the task objectives 
   were achieved, not on the number of failed attempts along the way.

Score each rubric category based on the agent's outputs. For EACH category, provide:
- score: the numeric score (0 to max points for that category)
- reasoning: a specific explanation of why this score was given, referencing the agent's outputs

Additionally, act as an EXECUTION JUDGE and report:
- execution_ok: true if the execution environment and tools were sufficiently functional to fairly judge the agent
- force_majeure: true if external factors (service outages, broken environment, etc.) substantially blocked the agent
- force_majeure_reason: a short explanation if force_majeure is true

Return JSON only:
{{{category_format}, "total": <0-{total_max}>, "execution_ok": <true/false>, "force_majeure": <true/false>, "force_majeure_reason": "<short explanation>"}}"""
    
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
    
    # Average each category - handles {"score": X, "reasoning": "..."} format
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
        category_scores[cat] = {
            "mean": statistics.mean(cat_scores),
            "scores": cat_scores,
            "reasoning": cat_reasoning
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


async def evaluate_output_async(task_dir: Path, output_dir: Path, client: AsyncOpenAI, 
                                 model: str, categories: dict, k: int) -> tuple[str, dict]:
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
    
    return output_dir.name, _aggregate_results(all_results, categories, failure_analysis)


async def evaluate_async(task: str, output_dirs: list[Path], k: int, eval_model: str, 
                         parallel_outputs: bool = False) -> tuple[dict, dict]:
    """Evaluate outputs with async parallelization.
    
    Args:
        task: Task name
        output_dirs: List of output directories to evaluate
        k: Number of evaluation rounds per output
        eval_model: Model to use for evaluation
        parallel_outputs: If True, evaluate multiple outputs in parallel.
                         If False, only parallelize the k rounds within each output.
    """
    client = AsyncOpenAI()
    task_dir = Path(__file__).parent / task
    _, categories = parse_rubric(task_dir / "rubric.txt")
    
    valid_dirs = [d for d in output_dirs if d.exists() and d.is_dir()]
    
    if parallel_outputs:
        # Parallel across all outputs AND all k rounds
        print(f"Running {len(valid_dirs)} outputs × {k} rounds in parallel...")
        tasks = [
            evaluate_output_async(task_dir, output_dir, client, eval_model, categories, k)
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
            
            results[output_dir.name] = _aggregate_results(all_results, categories, failure_analysis)
    
    return results, categories


def evaluate(task: str, output_dirs: list[Path], k: int, eval_model: str) -> dict:
    """Evaluate outputs k times (synchronous version for backwards compatibility)."""
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
        
        results[output_dir.name] = _aggregate_results(all_results, categories, failure_analysis)
    
    return results, categories


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
    parser.add_argument("--parallel", action="store_true", 
                        help="Parallelize across outputs (default: only parallelize k rounds per output)")
    parser.add_argument("--sync", action="store_true",
                        help="Use synchronous evaluation (no parallelization, for debugging)")
    
    args = parser.parse_args()
    task_dir = Path(__file__).parent / args.task
    agent_output = task_dir / "agent_output"
    
    # Get output directories
    if args.outputs:
        output_dirs = [agent_output / o if not Path(o).is_absolute() else Path(o) for o in args.outputs]
    elif agent_output.exists():
        output_dirs = sorted(agent_output.iterdir(), key=lambda p: p.name)
    else:
        print(f"No outputs found in {agent_output}")
        return
    
    print(f"Evaluating {len(output_dirs)} outputs ({args.k} rounds each)")
    
    if args.sync:
        print("Using synchronous evaluation (no parallelization)")
        results, categories = evaluate(args.task, output_dirs, args.k, args.model)
    else:
        mode = "fully parallel" if args.parallel else "parallel k rounds per output"
        print(f"Using async evaluation ({mode})")
        results, categories = asyncio.run(
            evaluate_async(args.task, output_dirs, args.k, args.model, args.parallel)
        )
    
    # Save
    with open(task_dir / f"{args.output}.json", "w") as f:
        json.dump({"categories": {k: v["max"] for k, v in categories.items()}, "results": results}, f, indent=2)
    
    plot_results(results, str(task_dir / f"{args.output}.png"))
    
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


if __name__ == "__main__":
    main()
