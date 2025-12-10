#!/usr/bin/env python3
"""
Evaluate agent outputs against a rubric.

Usage:
    python eval.py --task ph_sensitive_design --k 3
    python eval.py --task ph_sensitive_design --outputs run1 run2 --k 3
"""

import re
import json
import argparse
from pathlib import Path
import statistics
from dataclasses import dataclass, field
from typing import Optional

from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()
load_dotenv(Path(__file__).parent / ".env")


# ============================================================================
# Tool Failure Analysis
# ============================================================================

# Patterns indicating service unavailability (uncontrollable failures)
SERVICE_UNAVAILABLE_PATTERNS = [
    r"connection\s*(refused|reset|timed?\s*out)",
    r"(502|503|504)\s*(bad gateway|service unavailable|gateway timeout)",
    r"service\s*(unavailable|temporarily\s*unavailable)",
    r"api\s*(error|unavailable|rate\s*limit)",
    r"timeout\s*(error|exceeded)",
    r"network\s*(error|unreachable)",
    r"server\s*(error|unavailable|down)",
    r"quota\s*(exceeded|limit)",
    r"resource\s*exhausted",
    r"internal\s*server\s*error",
    r"ECONNREFUSED",
    r"ETIMEDOUT",
    r"ENOTFOUND",
]

# Patterns indicating tool usage/understanding errors
USAGE_ERROR_PATTERNS = [
    r"(invalid|missing|required)\s*(argument|parameter|field)",
    r"(type|validation)\s*error",
    r"(unexpected|unknown)\s*(argument|parameter|key)",
    r"ImportError|ModuleNotFoundError",
    r"(syntax|name|attribute)\s*error",
    r"KeyError|IndexError|TypeError",
    r"(file|path)\s*not\s*found",
    r"permission\s*denied",
    r"invalid\s*(input|format|syntax)",
]


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
    failed_calls: int = 0
    
    # Categorized failures
    service_failures: list = field(default_factory=list)  # Uncontrollable
    usage_failures: list = field(default_factory=list)    # Agent misunderstanding
    recovered_failures: list = field(default_factory=list)  # Failed then succeeded
    
    # Summary flags
    has_uncontrollable_failures: bool = False
    should_skip_eval: bool = False
    skip_reason: Optional[str] = None
    
    def to_dict(self) -> dict:
        return {
            "total_tool_calls": self.total_tool_calls,
            "successful_calls": self.successful_calls,
            "failed_calls": self.failed_calls,
            "service_failures": len(self.service_failures),
            "usage_failures": len(self.usage_failures),
            "recovered_failures": len(self.recovered_failures),
            "has_uncontrollable_failures": self.has_uncontrollable_failures,
            "should_skip_eval": self.should_skip_eval,
            "skip_reason": self.skip_reason,
        }


def classify_error(error_content: str) -> str:
    """Classify an error as 'service', 'usage', or 'unknown'."""
    error_lower = error_content.lower()
    
    # Check for service unavailability first (higher priority)
    for pattern in SERVICE_UNAVAILABLE_PATTERNS:
        if re.search(pattern, error_lower, re.IGNORECASE):
            return "service"
    
    # Check for usage/understanding errors
    for pattern in USAGE_ERROR_PATTERNS:
        if re.search(pattern, error_content, re.IGNORECASE):
            return "usage"
    
    return "unknown"


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
            tc.error_type = classify_error(result)
            analysis.failed_calls += 1
            
            if tc.error_type == "service":
                analysis.service_failures.append({
                    "tool": tc.tool_name,
                    "call_id": call_id,
                    "error_snippet": result[:500]
                })
            elif tc.error_type == "usage":
                analysis.usage_failures.append({
                    "tool": tc.tool_name,
                    "call_id": call_id,
                    "arguments": tc.arguments,
                    "error_snippet": result[:500]
                })
        else:
            analysis.successful_calls += 1
        
        # Track history for recovery detection
        if tc.tool_name not in tool_history:
            tool_history[tc.tool_name] = []
        tool_history[tc.tool_name].append((call_id, not is_error, tc.error_type))
    
    # Detect recoveries: tool failed (usage error) then succeeded
    for tool_name, history in tool_history.items():
        had_usage_failure = False
        for i, (call_id, success, error_type) in enumerate(history):
            if not success and error_type == "usage":
                had_usage_failure = True
            elif success and had_usage_failure:
                # Found a recovery!
                analysis.recovered_failures.append({
                    "tool": tool_name,
                    "description": f"Agent initially failed to use {tool_name} but eventually succeeded"
                })
                had_usage_failure = False  # Reset for next potential recovery
    
    # Determine if we should skip evaluation
    if analysis.service_failures:
        # Check if service failures are significant (>30% of calls or critical tools failed)
        service_failure_rate = len(analysis.service_failures) / max(analysis.total_tool_calls, 1)
        
        # Critical tools that if unavailable, should skip eval
        critical_tools = {"tamarind_submit_job", "tamarind_poll_results", "run_python", "esmfold", "proteinmpnn"}
        critical_service_failures = [
            f for f in analysis.service_failures 
            if any(ct in f["tool"].lower() for ct in critical_tools)
        ]
        
        if service_failure_rate > 0.3 or critical_service_failures:
            analysis.has_uncontrollable_failures = True
            analysis.should_skip_eval = True
            analysis.skip_reason = (
                f"Service unavailability detected: {len(analysis.service_failures)} service failures "
                f"out of {analysis.total_tool_calls} total calls. "
                f"Critical tools affected: {[f['tool'] for f in critical_service_failures]}"
            )
    
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
    """Build context string for the evaluator about tool failures."""
    if analysis.total_tool_calls == 0:
        return ""
    
    lines = []
    
    # Recovered failures - agent should NOT be penalized
    if analysis.recovered_failures:
        lines.append("RECOVERED ERRORS (do NOT penalize - agent learned and fixed these):")
        for rf in analysis.recovered_failures:
            lines.append(f"  - {rf['tool']}: {rf['description']}")
        lines.append("")
    
    # Service failures - uncontrollable
    if analysis.service_failures:
        lines.append("SERVICE UNAVAILABILITY (do NOT penalize - external service issues):")
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
    
    # Summary stats
    lines.append(f"TOOL CALL SUMMARY: {analysis.successful_calls}/{analysis.total_tool_calls} successful")
    
    return "\n".join(lines)


def evaluate_once(task_dir: Path, output_dir: Path, client: OpenAI, model: str, 
                  categories: dict, failure_analysis: Optional[ToolFailureAnalysis] = None) -> dict:
    """Single evaluation round."""
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
    
    # Build failure context if available
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
EVALUATION GUIDELINES FOR TOOL FAILURES:
1. SERVICE UNAVAILABILITY: If a tool/service was unavailable (API down, timeout, connection refused), 
   do NOT penalize the agent. These are external factors beyond the agent's control.
2. RECOVERED ERRORS: If the agent initially misunderstood a tool but eventually figured it out and 
   succeeded, do NOT penalize for the initial failures. Learning and recovery is expected behavior.
3. PERSISTENT MISUSE: If the agent consistently failed to use a tool correctly despite multiple 
   attempts (and never recovered), this MAY factor into the evaluation - but focus on whether the 
   agent achieved the task goals through alternative means.
4. FOCUS ON OUTCOMES: Primarily evaluate based on the final outputs and whether the task objectives 
   were achieved, not on the number of failed attempts along the way.

Score each rubric category based on the agent's outputs. For EACH category, provide:
- score: the numeric score (0 to max points for that category)
- reasoning: a specific explanation of why this score was given, referencing the agent's outputs

Return JSON only:
{{{category_format}, "total": <0-{total_max}>}}"""

    response = client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        response_format={"type": "json_object"}
    )
    
    try:
        return json.loads(response.choices[0].message.content)
    except:
        return {"total": 0, "error": "Parse failed"}


def evaluate(task: str, output_dirs: list[Path], k: int, eval_model: str) -> dict:
    """Evaluate outputs k times."""
    client = OpenAI()
    task_dir = Path(__file__).parent / task
    _, categories = parse_rubric(task_dir / "rubric.txt")
    results = {}
    
    for output_dir in output_dirs:
        if not output_dir.exists():
            print(f"Skip: {output_dir} not found")
            continue
        
        # Skip non-directories (e.g., .DS_Store)
        if not output_dir.is_dir():
            print(f"Skip: {output_dir} is not a directory")
            continue
        
        print(f"\n{output_dir.name}")
        
        # Analyze tool failures from agent log
        agent_log_path = output_dir / "agent_log.json"
        failure_analysis = analyze_tool_failures(agent_log_path)
        
        # Print failure analysis summary
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
        
        # Check if we should skip evaluation due to service unavailability
        if failure_analysis.should_skip_eval:
            print(f"  ⛔ SKIPPING EVALUATION: {failure_analysis.skip_reason}")
            results[output_dir.name] = {
                "scores": [],
                "mean": None,
                "std": None,
                "breakdown": {},
                "reasoning": [],
                "skipped": True,
                "skip_reason": failure_analysis.skip_reason,
                "failure_analysis": failure_analysis.to_dict(),
            }
            continue
        
        scores = []
        all_results = []
        
        for i in range(k):
            result = evaluate_once(task_dir, output_dir, client, eval_model, categories, failure_analysis)
            all_results.append(result)
            score = result.get("total", 0)
            scores.append(score)
            print(f"  Round {i+1}: {score}")
        
        # Average each category - now handles {"score": X, "reasoning": "..."} format
        category_scores = {}
        for cat in categories:
            # Extract scores and reasoning from each round
            cat_data = [r.get(cat, {}) for r in all_results]
            # Handle both old format (just number) and new format ({"score": X, "reasoning": "..."})
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
                "reasoning": cat_reasoning  # Per-category reasoning for each round
            }
        
        results[output_dir.name] = {
            "scores": scores,
            "mean": statistics.mean(scores),
            "std": statistics.stdev(scores) if len(scores) > 1 else 0,
            "breakdown": category_scores,
            "skipped": False,
            "failure_analysis": failure_analysis.to_dict(),
        }
    
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
    
    results, categories = evaluate(args.task, output_dirs, args.k, args.model)
    
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
