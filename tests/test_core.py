"""Essential tests for core modules (no API calls)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "tasks"))

import pytest

# === Import Tests ===
def test_imports():
    """All core modules import without errors."""
    from tamarind_client import TamarindClient, _load_env_file
    from dag_tracker import DAGTracker, Node, Edge
    from prompts import build_agent_system_prompt, build_evaluator_prompt
    from eval import parse_rubric, ToolFailureAnalysis

# === DAGTracker Tests ===
def test_dag_tracker_nodes_edges():
    """DAGTracker creates nodes and edges correctly."""
    from dag_tracker import DAGTracker
    dag = DAGTracker()
    a1 = dag.add_artifact("input.pdb")
    f1 = dag.add_function("process")
    a2 = dag.add_artifact("output.pdb")
    dag.add_edge(a1, f1)
    dag.add_edge(f1, a2)
    assert len(dag.nodes) == 3
    assert len(dag.edges) == 2

def test_dag_tracker_serialization(tmp_path):
    """DAGTracker save/load roundtrip."""
    from dag_tracker import DAGTracker
    dag = DAGTracker()
    dag.add_artifact("test.pdb")
    path = tmp_path / "dag.json"
    dag.save(path)
    loaded = DAGTracker.load(path)
    assert len(loaded.nodes) == 1

# === Rubric Parsing Tests ===
def test_parse_rubric(tmp_path):
    """parse_rubric extracts categories and points."""
    from eval import parse_rubric
    rubric = tmp_path / "rubric.txt"
    rubric.write_text("Criterion 1: Analysis\nMax Points: 25\n\nCriterion 2: Results\nMax Points: 50\n")
    text, cats = parse_rubric(rubric)
    assert "analysis" in cats
    assert cats["analysis"]["max"] == 25
    assert cats["results"]["max"] == 50

# === Prompt Building Tests ===
def test_build_agent_prompt():
    """Agent prompt includes task description."""
    from prompts import build_agent_system_prompt
    prompt = build_agent_system_prompt("Design a protein binder")
    assert "Design a protein binder" in prompt
    assert "AVAILABLE TOOLS" in prompt

def test_build_evaluator_prompt():
    """Evaluator prompt includes all sections."""
    from prompts import build_evaluator_prompt
    prompt = build_evaluator_prompt(
        question="Test task", rubric_text="Test rubric", output_dir_name="run1",
        files_json="[]", contents_json="{}", failure_context="",
        category_format='"test": {"score": 10}', total_max=100
    )
    assert "Test task" in prompt
    assert "Test rubric" in prompt
    assert "execution_ok" in prompt

# === ToolFailureAnalysis Tests ===
def test_tool_failure_analysis():
    """ToolFailureAnalysis serializes correctly."""
    from eval import ToolFailureAnalysis
    analysis = ToolFailureAnalysis(total_tool_calls=5, successful_calls=3, failed_calls=2)
    d = analysis.to_dict()
    assert d["total_tool_calls"] == 5
    assert d["successful_calls"] == 3

