# pH-Sensitive Protein Design Task

## Background

Histidine's imidazole side chain has a pKa of approximately 6.0, which falls between physiological pH (7.4) and endosomal/tumor microenvironment pH (5.5-6.0). This makes Histidine ideal for pH-sensing applications:

- At neutral pH (7.4): His is deprotonated (neutral) and can form hydrogen bonds
- At acidic pH (~6.0): His becomes protonated (positive), disrupting hydrogen bonds

For pH sensing to work effectively, the His network must be **buried in the protein core**. Surface-exposed His residues have normal pKa (~6.0), but buried His residues experience pKa shifts due to the hydrophobic environment. When buried His residues are protonated at low pH, this disrupts core hydrogen bonds and destabilizes the protein fold.

## Task

Using PDB ID **5L33** (a de novo designed NTF2 fold) as your scaffold, implement the following workflow:

### Step 1: Identify Buried Core Residues
Calculate Solvent Accessible Surface Area (SASA) for each residue using BioPython's `Bio.PDB.SASA.ShrakeRupley`. Compute relative SASA by dividing by the maximum SASA for each amino acid type. Select residues with **relative SASA < 0.25** as core residues.

Reference MAX_SASA values (Å²):
```
ALA: 129, ARG: 274, ASN: 195, ASP: 193, CYS: 167, GLN: 225, GLU: 223,
GLY: 104, HIS: 224, ILE: 197, LEU: 201, LYS: 236, MET: 224, PHE: 240,
PRO: 159, SER: 155, THR: 172, TRP: 285, TYR: 263, VAL: 174
```

### Step 2: Find Network Positions
Calculate Cβ-Cβ distances between core residue pairs (use Cα for Glycine). Select pairs suitable for His-containing hydrogen bond networks:
- His-His pairs: Cβ-Cβ distance **5.5-8.5 Å**
- His-Arg or His-Lys pairs: Cβ-Cβ distance **6.0-10.0 Å**

### Step 3: Design Sequences with ProteinMPNN
Use ProteinMPNN via Tamarind API to redesign the sequence:
- Upload the scaffold PDB
- Constrain network positions to His (or Arg/Lys) using `bias_AA_per_residue`
- Generate at least 1 sequence

### Step 4: Validate with Structure Prediction
Use ESMFold via Tamarind API:
- Predict structure for at least 1 designed sequence
- Extract pLDDT scores from B-factor column
- Filter predictions with **mean pLDDT > 80**
- **Verify Motif**: Confirm that the predicted structure contains the intended motif residues (e.g., His) at the correct network positions.

## Deliverables

- Working Python script implementing the workflow
- Output JSON files with core residues, network positions, and designed sequences
- Predicted structure files

## Available Tools

- **BioPython**: PDB parsing, SASA calculation
- **Tamarind API**: ProteinMPNN, ESMFold
- **NumPy**: Distance calculations
