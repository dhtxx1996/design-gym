# pH-Sensitive Protein Design Task

## Problem Statement

You are tasked with engineering a pH-sensitive "switch" into a stable protein scaffold. The goal is to introduce a buried Histidine-mediated hydrogen bond network that is stable at neutral pH (7.4) but becomes destabilizing at acidic pH (~6.0), triggering a conformational change or unfolding.

Using the provided scaffold `data/scaffold.pdb` (a de novo designed NTF2 fold, PDB ID 5L33), develop a computational workflow to:

1.  **Identify Core Residues**: Analyze the structure to find buried residues suitable for core engineering. pH-sensing residues must be buried to experience the pKa shift necessary for sensing physiological pH changes.
2.  **Select Network Positions**: Identify a geometric motif (e.g., a pair or triplet of residues) within the core that can be mutated to Histidines (or His-Polar pairs) to form a hydrogen bond network.
3.  **Computational Design**:
    *   Mutate the selected positions to Histidine.
    *   Redesign the surrounding sequence to pack against the new Histidines and stabilize the neutral-pH state (Inverse Folding).
4.  **Validation**:
    *   Predict the structure of your designed sequences (Forward Folding).
    *   Verify that the predicted structure adopts the target fold (high pLDDT/low RMSD).
    *   **Crucially**, verify that the intended Histidine network is formed in the predicted structure (check geometry/distances).

## Deliverables

1.  **Data Files** (saved to output directory):
    *   `core.json`: List of identified buried residues with their properties.
    *   `network.json`: The selected positions for the histidine network with geometric data.
    *   `designs.json`: The designed sequences.
    *   `predictions.json`: Validation metrics (pLDDT, etc.) for the designs.
2.  **Analysis**: A brief explanation of why you chose the specific positions.

## Technical Hints

### SASA Calculation
- Use BioPython's `ShrakeRupley` class (NOT DSSP which requires external tools):
  ```python
  from Bio.PDB.SASA import ShrakeRupley
  sr = ShrakeRupley()
  sr.compute(structure, level="R")  # Per-residue SASA
  # Access via: residue.sasa
  ```
- Relative SASA = absolute SASA / max SASA for residue type
- Max SASA values: ALA~129, LEU~201, PHE~240, etc. (look up standard values)
- Buried residues: relative SASA < 0.25

### PDB Numbering
- **IMPORTANT**: The scaffold PDB has non-standard residue numbering (starts at -1, 0, 1...)
- Always use actual PDB residue IDs from `residue.id[1]`, not Python list indices
- When using Tamarind tools, residue identifiers must be in format: `{chain}{residue_number}` (e.g., "A5", "A-1")

### Geometry for His Networks
- Optimal Cβ-Cβ distance for His-His hydrogen bonds: 5.5 - 8.0 Å (optimal ~6.5 Å)
- Use Cα for Glycine (no Cβ)

### ProteinMPNN via Tamarind
- First upload PDB with `tamarind_upload_file`
- Use `tamarind_get_tool_spec("proteinmpnn")` to see all parameters
- Key parameters:
  - `pdbFile`: filename of uploaded PDB
  - `numSequences`: number of designs to generate
  - `temperature`: diversity (0.1 = low diversity, 0.3 = higher)
  - `bias_AA_per_residue`: JSON string to bias specific amino acids at positions
    - Format: `'{"A5": {"H": 10.0}, "A9": {"H": 10.0}}'`
    - High positive values (e.g., 10.0) strongly bias toward that amino acid
- **Note**: Do NOT use `designedChains` or `verifySequences` together with `bias_AA_per_residue` (API limitations)

### ESMFold via Tamarind
- Use `tamarind_get_tool_spec("esmfold")` to see parameters
- Submit sequence directly: `{"sequence": "MKTVRQ..."}`
- Output PDB has pLDDT scores in the B-factor column

## Grading Criteria

Your solution will be evaluated on:
*   **Conceptual Understanding**: Correct identification of buried residues and appropriate motif selection.
*   **Geometric Rigor**: Evidence that the installed motif satisfies geometric constraints for hydrogen bonding.
*   **Workflow Logic**: Proper use of inverse folding (design) followed by forward folding (validation).
*   **Success Metrics**: Generation of a design that is predicted to fold into the target structure with the installed network intact.
