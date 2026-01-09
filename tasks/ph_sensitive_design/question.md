# pH-Sensitive Protein Design Task

## Problem Statement

Engineer a pH-sensitive "switch" into a stable protein scaffold. The goal is to introduce a buried Histidine-mediated hydrogen bond network that is stable at neutral pH (7.4) but becomes destabilizing at acidic pH (~6.0), triggering a conformational change or unfolding.

Using PDB ID **5L33** (a de novo designed NTF2 fold) as your scaffold, develop a computational workflow that:

1. **Identifies buried core residues** suitable for installing pH-sensing networks
2. **Discovers positions** where His-containing hydrogen bond networks can be installed
3. **Designs sequences** with the network installed and surrounding residues repacked
4. **Validates designs** through structure prediction, confirming the network forms correctly

## Success Criteria

A successful solution will:
- Produce designed sequences containing buried His networks
- Demonstrate high-confidence structure predictions (pLDDT > 80)
- Show geometric evidence that the His network is formed in predicted structures
- Explain conceptually why Histidine is chosen for pH sensing
- Explain why the His network must be buried in the protein core
- Discuss potential biological applications of such pH-sensitive proteins

## Background

- Histidine's pKa (~6.0) makes it ideal for physiological pH sensing
- Buried residues experience pKa shifts that enable pH sensing
- At low pH, His protonation disrupts hydrogen bonds, destabilizing the core

## Deliverables

- Working Python script implementing the pipeline
- Output files with core residue analysis, network positions, designed sequences, and structure predictions
- Written explanations of the design rationale and biochemistry
