# pH-Sensitive Protein Design Task

Design a pH-sensitive protein by installing buried Histidine networks that destabilize at acidic pH.

## Files

| File | Description |
|------|-------------|
| `question_{easy,medium,hard}.md` | Task prompts at different difficulty levels |
| `workflow.py` | Reference implementation (100/100 on rubric) |
| `rubric.txt` | Scoring criteria (8 criteria, 100 points) |
| `data/scaffold.pdb` | Input structure (PDB 5L33) |

## Question Difficulty Levels

| Level | Guidance Provided | Agent Must Infer |
|-------|-------------------|------------------|
| **Easy** | Step-by-step workflow, exact parameters (SASA<0.25, Cb-Cb 5.5-8.5Å), tool names, reference values | Minimal |
| **Medium** | High-level workflow steps, success criteria, hints, background concepts | Parameter values, tool selection |
| **Hard** | Goal and deliverables only | Entire approach, parameters, tools, validation strategy |

## Quick Start

```bash
python workflow.py --pdb data/scaffold.pdb --output output/
```
