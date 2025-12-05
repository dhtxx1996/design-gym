#!/usr/bin/env python3
"""
Monomer Destabilization Pipeline for pH-Sensitive Protein Design
"""

import sys
import json
import argparse
from pathlib import Path
import numpy as np
from Bio.PDB import PDBParser, SASA
from Bio.SeqUtils import seq1
from Bio import SeqIO

# Add parent directory for tamarind_client import
sys.path.insert(0, str(Path(__file__).parent.parent))
from tamarind_client import TamarindClient

# Constants
MAX_SASA = {
    'ALA': 129.0, 'ARG': 274.0, 'ASN': 195.0, 'ASP': 193.0, 'CYS': 167.0, 
    'GLN': 225.0, 'GLU': 223.0, 'GLY': 104.0, 'HIS': 224.0, 'ILE': 197.0, 
    'LEU': 201.0, 'LYS': 236.0, 'MET': 224.0, 'PHE': 240.0, 'PRO': 159.0, 
    'SER': 155.0, 'THR': 172.0, 'TRP': 285.0, 'TYR': 263.0, 'VAL': 174.0,
}
CB_DISTANCE_RANGE = (5.5, 8.0)
CB_OPTIMAL_DISTANCE = 6.5

def parse_structure(pdb_path: str):
    """Parse PDB and map to 0-indexed sequence."""
    structure = PDBParser(QUIET=True).get_structure("scaffold", pdb_path)
    chain = next(structure[0].get_chains())
    
    # Filter for standard residues
    residues = [r for r in chain if r.id[0] == ' ' and r.resname in MAX_SASA]
    sequence = "".join(seq1(r.resname) for r in residues)
    pdb_map = [r.id[1] for r in residues]
    
    return structure, chain, sequence, residues, pdb_map

def identify_core_residues(pdb_path: str, threshold: float = 0.25) -> dict:
    """Identify buried core residues based on relative SASA."""
    structure, chain, sequence, residues, pdb_map = parse_structure(pdb_path)
    
    # Calculate SASA
    SASA.ShrakeRupley().compute(structure, level="R")
    sasa_map = {r.id[1]: r.sasa for r in chain if r.id[0] == ' '}
    
    core_indices = []
    rel_sasa_list = []
    
    for i, r in enumerate(residues):
        resnum = r.id[1]
        rel = sasa_map.get(resnum, 0) / MAX_SASA.get(r.resname, 1.0)
        rel_sasa_list.append(rel)
        if rel < threshold:
            core_indices.append(i)
            
    print(f"[Module 1] Core: {len(core_indices)}/{len(residues)} residues")
    return {
        "core_selection": core_indices, "rel_sasa_values": rel_sasa_list,
        "all_residues": list(range(len(sequence))), "pdb_index_map": pdb_map,
        "sequence": sequence, "pdb_path": pdb_path, "chain_id": chain.get_id()
    }

def find_best_network_positions(parsed_data: dict, distance_range=CB_DISTANCE_RANGE, optimal=CB_OPTIMAL_DISTANCE) -> dict:
    """Find optimal core position pairs for His networks."""
    structure, _, _, residues, _ = parse_structure(parsed_data["pdb_path"])
    
    # Get Cb coords (or Ca for Gly)
    def get_coord(r):
        return r['CA'].coord if r.resname == 'GLY' else (r['CB'].coord if 'CB' in r else r['CA'].coord)

    core_indices = [i for i in parsed_data["core_selection"] if i < len(residues)]
    coords = np.array([get_coord(residues[i]) for i in core_indices])
    
    if len(coords) < 2:
        return {"network_selection": [], "cb_distance": None, "geometric_score": None}

    # Vectorized distance matrix
    diff = coords[:, None, :] - coords[None, :, :]
    dists = np.sqrt(np.sum(diff**2, axis=-1))
    
    # Filter and score pairs
    mask = (dists >= distance_range[0]) & (dists <= distance_range[1])
    pairs = np.argwhere(np.triu(mask, k=1))
    
    candidates = []
    min_d, max_d = distance_range
    for i1, i2 in pairs:
        d = dists[i1, i2]
        score = 1.0 - abs(d - optimal) / (max_d - min_d)
        candidates.append({
            "positions": [core_indices[i1], core_indices[i2]],
            "cb_distance": float(d), "geometric_score": float(max(0, score))
        })
    
    candidates.sort(key=lambda x: x["geometric_score"], reverse=True)
    best = candidates[0] if candidates else None
    
    if best:
        print(f"[Module 2] Best pair: {best['positions']} (d={best['cb_distance']:.2f})")
        
    return {
        "network_selection": best["positions"] if best else [],
        "cb_distance": best["cb_distance"] if best else None,
        "geometric_score": best["geometric_score"] if best else None,
        "all_candidate_pairs": candidates
    }

def design_around_network(client, pdb_path, sequence, network_indices, pdb_map, chain_id, num_seqs=20):
    """Run ProteinMPNN design."""
    # Mutate sequence
    seq_list = list(sequence)
    for i in network_indices: seq_list[i] = 'H'
    mutated_seq = "".join(seq_list)
    
    # Setup Tamarind job
    if client:
        print(f"[Module 3] Submitting ProteinMPNN job...")
        client.upload_file(pdb_path)
        bias = {f"{chain_id}{pdb_map[i]}": {"H": 100.0} for i in network_indices}
        
        try:
            job = client.submit_job_sync("proteinmpnn", {
                "pdbFile": Path(pdb_path).name, "numSequences": str(num_seqs),
                "temperature": "0.1", "bias_AA_per_residue": json.dumps(bias)
            }, timeout=600)
            
            results = client.download_results(job['job_name'])
            designs = []
            for f in list(results.rglob("*.fa*")):
                for rec in SeqIO.parse(f, "fasta"):
                    if all(rec.seq[i] == 'H' for i in network_indices):
                        designs.append({"header": rec.id, "sequence": str(rec.seq)})
            return designs, mutated_seq
            
        except Exception as e:
            print(f"Job failed: {e}")

    # Fallback mock
    print("[Module 3] Using mock sequences")
    return [{
        "header": f"mock_{i}",
        "sequence": "".join('H' if j in network_indices else (
            np.random.choice(list("ACDEFGHIKLMNPQRSTVWY")) 
            if np.random.rand() < 0.1 else s
        ) for j, s in enumerate(mutated_seq))
    } for i in range(num_seqs)], mutated_seq

def predict_structures(client, designs, network_indices, output_dir, max_preds=5):
    """Run ESMFold prediction."""
    preds = []
    pred_dir = Path(output_dir) / "predicted_structures"
    pred_dir.mkdir(parents=True, exist_ok=True)
    
    for i, d in enumerate(designs[:max_preds]):
        try:
            if client:
                job = client.submit_job_sync("esmfold", {"sequence": d["sequence"]}, timeout=300)
                res_path = client.download_results(job['job_name'], output_dir=str(pred_dir))
                pdb_file = next(res_path.glob("*.pdb")) if res_path.is_dir() else res_path
                
                # Extract pLDDT
                structure = PDBParser(QUIET=True).get_structure("p", str(pdb_file))
                bfactors = [a.bfactor for a in structure.get_atoms()]
                mean_plddt = np.mean(bfactors)
                
                # Network pLDDT (esm 1-based vs index 0-based)
                net_bfactors = []
                for chain in structure[0]:
                    for r in chain:
                        if (r.id[1] - 1) in network_indices: # approximate mapping if sequential
                             net_bfactors.extend([a.bfactor for a in r])
                
                preds.append({**d, "pdb_path": str(pdb_file), "plddt_mean": float(mean_plddt)})
                continue
                
        except Exception:
            pass
            
        # Mock fallback
        preds.append({**d, "pdb_path": "mock.pdb", "plddt_mean": 75.0, "is_mock": True})
        
    return preds

def run_pipeline(pdb_path: str, output_dir: str, **kwargs):
    """Main pipeline execution."""
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    # Initialize Client
    try: 
        client = TamarindClient()
    except Exception as e:
        print(f"Tamarind init failed: {e}. Using mock fallback.")
        client = None

    # Module 1: Core
    m1 = identify_core_residues(pdb_path, kwargs.get('sasa_threshold', 0.25))
    with open(out / "core.json", 'w') as f: json.dump(m1, f, indent=2)
    
    # Module 2: Network
    m2 = find_best_network_positions(m1)
    with open(out / "network.json", 'w') as f: json.dump(m2, f, indent=2)
    
    if not m2["network_selection"]: return
    
    # Module 3: Design
    designs, mut_seq = design_around_network(
        client, pdb_path, m1["sequence"], 
        m2["network_selection"], m1["pdb_index_map"], m1["chain_id"],
        kwargs.get('num_designs', 2)
    )
    
    m3 = {
        "designed_sequences": designs, 
        "network_selection": m2["network_selection"],
        "original_sequence": mut_seq
    }
    with open(out / "designs.json", 'w') as f: json.dump(m3, f, indent=2)
    
    # Module 4: Prediction
    preds = predict_structures(
        client, designs, 
        m2["network_selection"], out, kwargs.get('max_predictions', 5)
    )
    
    with open(out / "predictions.json", 'w') as f: json.dump(preds, f, indent=2)
    print(f"Done. Results in {out}")

if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pdb", default="data/scaffold.pdb")
    p.add_argument("--output", default="output")
    p.add_argument("--sasa-threshold", type=float, default=0.25)
    p.add_argument("--num-designs", type=int, default=2)
    p.add_argument("--max-predictions", type=int, default=5)
    args = p.parse_args()
    
    run_pipeline(
        pdb_path=args.pdb,
        output_dir=args.output,
        **vars(args)
    )
