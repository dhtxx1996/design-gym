- how to determine interface?
- find an example binder to make the scoring redesign easy
- can we use a pyRosetta alternative 
- use the proteinMPNN in the Biomni sandbox?
- impelement custom Rosetta HBNet? (but it requires Rotamer sampling tho, https://docs.rosettacommons.org/docs/latest/scripting_documentation/RosettaScripts/Movers/movers_pages/HBNetMover)


- Given a few candidate structures and determine which has the PH senstivie properties 

- Key point: The benchmark should not be dependent on a specific tool or environment; even listing tool options or specifying tools may bias the agent's choices and compromise the neutrality of the evaluation.

- potential PDBs with PH sensitive interface
	•	5DD8 – HucR histidine stack at dimer interface (natural pH switch).  ￼
	•	6KMC – engineered Protein G mutant with His-mediated pH-dependent IgG binding.  ￼
	•	5WLE – PPS PHD finger, His in aromatic cage controlling pH-dependent H3K4me3 binding.  ￼
	•	5GU6 / 5XWM – ERp44, His cluster + Zn²⁺ regulating pH-dependent dimerization/client binding.  ￼
	•	8YG6 / 8YG8 / 3DUZ – GP64 viral fusion protein, multi-His pH-sensitive trimer interface.


