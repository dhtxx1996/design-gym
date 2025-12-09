# pH-Sensitive Protein Design Task

## Problem Statement

Engineer a pH-sensitive "switch" into a stable protein scaffold. The goal is to introduce a buried Histidine-mediated hydrogen bond network that is stable at neutral pH (7.4) but becomes destabilizing at acidic pH (~6.0), triggering a conformational change or unfolding.

Use the de novo designed NTF2 fold **PDB ID 5L33** as your scaffold. Download it to your working directory:
```python
urllib.request.urlretrieve("https://files.rcsb.org/download/5L33.pdb", "scaffold.pdb")
```

Develop a computational workflow that:

1. **Identifies buried core residues** suitable for installing pH-sensing networks
2. **Discovers positions** where His-containing hydrogen bond networks can be installed
3. **Designs sequences** with the network installed and surrounding residues repacked
4. **Validates designs** through structure prediction, confirming the network forms correctly

## Success Criteria

A successful solution will:
- Produce designed sequences containing buried His networks
- Demonstrate high-confidence structure predictions (pLDDT > 80)
- Show geometric evidence that the His network is formed in predicted structures
- **Explain conceptually why Histidine is chosen** for pH sensing (hint: consider its pKa ~6.0 relative to physiological/endosomal pH ranges)
- **Explain why the His network must be buried** in the protein core (hint: consider pKa shifts in buried vs. surface environments)
- **Discuss potential biological applications** of such pH-sensitive proteins (e.g., endosomal release for drug delivery, biosensors, pH-triggered conformational switches)

## Background

- Histidine's pKa (~6.0) makes it ideal for physiological pH sensing
- Buried residues experience pKa shifts that enable pH sensing
- At low pH, His protonation disrupts hydrogen bonds, destabilizing the core
