"""Minimal tests for TamarindClient API."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent / "tasks"))

import pytest
from tamarind_client import TamarindClient

@pytest.fixture
def client():
    return TamarindClient()

def test_list_tools(client):
    """API connectivity - list available tools."""
    tools = client.list_tool_names()
    assert "esmfold" in tools
    assert "proteinmpnn" in tools

def test_get_tool_spec(client):
    """Tool discovery - get ESMFold spec."""
    spec = client.get_tool_spec("esmfold")
    assert spec is not None
    assert "sequence" in str(spec).lower()

def test_upload_file(client, tmp_path):
    """File upload."""
    pdb = tmp_path / "test.pdb"
    pdb.write_text("ATOM      1  CA  ALA A   1       0.0   0.0   0.0  1.00  0.00           C\nEND\n")
    result = client.upload_file(str(pdb))
    assert "filename" in result or "error" not in str(result).lower()

# @pytest.mark.slow
# def test_esmfold_job(client):
#     """Submit ESMFold job, poll, verify completion."""
#     result = client.submit_job_sync("esmfold", {"sequence": "GIVEQCCTSICSLYQLENYCN"}, timeout=300)
#     assert result.get("job_name")
#     assert result.get("final_status", {}).get("status") in ["completed", "pending"]

# @pytest.mark.slow  
# def test_proteinmpnn_job(client):
#     """Submit ProteinMPNN job with bias."""
#     import json
#     # Use a known file on Tamarind
#     params = {
#         "pdbFile": "9f6o.pdb",
#         "numSequences": "2",
#         "bias_AA_per_residue": json.dumps({"A1": {"G": 1.0}})
#     }
#     result = client.submit_job_sync("proteinmpnn", params, timeout=300)
#     assert result.get("job_name")
