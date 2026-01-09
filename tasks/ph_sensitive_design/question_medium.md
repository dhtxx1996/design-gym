# pH-Sensitive Protein Design Task

## Problem Statement

Engineer a pH-sensitive "switch" into a stable protein scaffold. The goal is to introduce a buried Histidine-mediated hydrogen bond network that is stable at neutral pH (7.4) but becomes destabilizing at acidic pH (~6.0), triggering a conformational change or unfolding.

Using PDB ID **5L33** (a de novo designed NTF2 fold) as your scaffold, develop a computational workflow that:

1. **Identifies buried core residues** suitable for installing pH-sensing networks
2. **Discovers positions** where His-containing hydrogen bond networks can be installed
3. **Designs sequences** with the network installed and surrounding residues repacked
4. **Validates designs** through structure prediction, confirming the network forms correctly

## Scientific Requirements

Your solution must demonstrate understanding of the following:

1. **Why Histidine for pH sensing?**
   Explain why Histidine is the amino acid of choice for this application. Consider its pKa relative to the target pH range.

2. **Why must the network be buried?**
   Explain the relationship between residue burial and pKa shifts. Why wouldn't a surface-exposed His network work for this application?

3. **Network geometry considerations**
   Justify your criteria for selecting network positions. What structural features make a good His hydrogen bond network?

4. **Biological applications**
   Discuss potential applications of pH-sensitive proteins (e.g., drug delivery, biosensing).

## Technical Hints

- Buried residues typically have relative SASA < 0.25
- His-His hydrogen bonds require appropriate Cβ-Cβ geometry (5-9 Å range)
- Consider His-Arg and His-Lys networks as alternatives to His-His
- Ensure network residues are separated in sequence (e.g., >10 residues apart) to form tertiary contacts
- Verify geometric orientation: Ca-Cb vectors should point towards each other (positive dot product with connecting vector)
- Use inverse folding tools to redesign around fixed network positions
- Filter structure predictions by confidence scores (pLDDT) and RMSD to the scaffold
- Explicitly verify that the designed motif residues (identity and location) are preserved in the predicted structure

## Deliverables

- Working Python script implementing the complete workflow
- Output files demonstrating results at each step
- Written explanations addressing the scientific requirements above

## Available Tools

- **BioPython**: PDB parsing, SASA calculation
- **Tamarind API**: ProteinMPNN, ESMFold, AlphaFold, and other computational biology tools
- **NumPy**: Numerical calculations
