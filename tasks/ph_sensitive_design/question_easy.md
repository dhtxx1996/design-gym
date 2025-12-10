# pH-Sensitive Protein Design Task (EASY)

## Problem Statement

Engineer a pH-sensitive "switch" into a stable protein scaffold. The goal is to introduce a buried Histidine-mediated hydrogen bond network that is stable at neutral pH (7.4) but becomes destabilizing at acidic pH (~6.0), triggering a conformational change or unfolding.

Use the de novo designed NTF2 fold **PDB ID 5L33** as your scaffold. Download it to your working directory:
```python
urllib.request.urlretrieve("https://files.rcsb.org/download/5L33.pdb", "scaffold.pdb")
```

Implement the following workflow:

---

## Step-by-Step Workflow

### Step 1: Identify Buried Core Residues
**Goal**: Find residues that are buried inside the protein (not surface-exposed).

**Method**: Calculate Solvent Accessible Surface Area (SASA) for each residue.
- Use BioPython's `Bio.PDB.SASA.ShrakeRupley` to compute SASA
- Calculate **relative SASA** = (residue SASA) / (max SASA for that amino acid type)
- Select residues with **relative SASA < 0.25** as "core" residues

**Reference MAX_SASA values** (Å²):
```
ALA: 129, ARG: 274, ASN: 195, ASP: 193, CYS: 167, GLN: 225, GLU: 223,
GLY: 104, HIS: 224, ILE: 197, LEU: 201, LYS: 236, MET: 224, PHE: 240,
PRO: 159, SER: 155, THR: 172, TRP: 285, TYR: 263, VAL: 174
```

### Step 2: Find His Network Positions
**Goal**: Identify pairs of core residues suitable for His-His (or His-Arg, His-Lys) hydrogen bond networks.

**Method**: Calculate Cβ-Cβ distances between core residue pairs.
- For Glycine (no Cβ), use Cα instead
- Select pairs with Cβ-Cβ distance in **5.5 - 8.5 Å** range for His-His networks
- For His-Arg or His-Lys pairs, use **6.0 - 10.0 Å** range
- Choose the pair closest to optimal distance (~6.5 Å for His-His)

### Step 3: Design Sequences with ProteinMPNN
**Goal**: Redesign the protein sequence while constraining the His network positions.

**Method**: Use ProteinMPNN via Tamarind API:
- Upload the scaffold PDB
- Constrain network residues using `bias_AA_per_residue` with positive bias (e.g., 5.0-100.0) to favor His (or Arg/Lys) at network positions
- Generate **at least 1 sequence** with temperature = 0.1
- Verify output sequences have correct residues at network positions

### Step 4: Validate with Structure Prediction
**Goal**: Confirm designed sequences fold correctly with the His network formed.

**Method**: Use ESMFold via Tamarind API:
- Predict structures for **at least 1 design**
- Extract pLDDT scores from B-factor column of output PDB
- **Filter**: Keep only predictions with **mean pLDDT > 80**
- Verify His network geometry in predicted structures

---

## Required Explanations

Your solution must include explanations for:

1. **Why Histidine?**
   - Histidine's imidazole side chain has pKa ~6.0
   - This falls between physiological pH (7.4) and endosomal pH (5.5-6.0)
   - At neutral pH: His is deprotonated (neutral) → forms H-bonds
   - At acidic pH: His is protonated (positive) → H-bonds disrupted

2. **Why must the network be buried?**
   - Surface His have normal pKa (~6.0)
   - Buried His experience pKa shifts due to hydrophobic environment
   - Protonation of BURIED His disrupts core H-bonds → destabilizes protein
   - Surface His protonation has minimal structural effect

3. **Biological Applications**:
   - Endosomal release for drug delivery
   - Tumor-targeting therapeutics
   - pH biosensors
   - Controlled protein switches

---

## Available Tools

- **BioPython**: PDB parsing, SASA calculation, PDB download
- **Tamarind API**: ProteinMPNN (inverse folding), ESMFold (structure prediction), and other tools
- **NumPy**: Distance calculations

