# Methods · Interpreting structural evidence

[← Back to the project](../README.md)

This guide explains the structural metrics and the parameters to retain when comparing models. The analyses operate on the input structural models and support the development and prioritization of research hypotheses.

## Cross-region contacts

For residue sets A and B in the selected chain, the program calculates the minimum Euclidean distance between atoms in each residue pair. Atom selection follows `_get_residue_atom_coords` in `PP.py`; the current implementation should not be described as a heavy-atom-only calculation.

| Metric | Definition and interpretation |
| :--- | :--- |
| `CrossContactPairs` | Number of residue pairs whose minimum atomic distance is within the selected cutoff |
| `CrossContactDensity` | Contact count divided by the theoretical residue-pair count, `len(A) × len(B)` |
| `CrossContactMinDist` | Smallest observed atomic distance between the two sets |
| `Pairs@…` | Contact counts at different distance thresholds |
| Top-K distances | Distances and averages for the closest residue pairs, complementing integer counts |
| Differences from WT | Distance changes for matching residue pairs relative to the baseline model |

The cutoff, set size, and input coordinates affect the results. Counts are coarse for small sets; strict thresholds can produce all-zero counts, while permissive thresholds can saturate. Report minimum distances, a documented set of cutoffs, and the relevant structural figures together.

The program can also select a cutoff based on discrimination between models. Treat that choice as exploratory: describe the selection method and show the threshold scan so readers can assess its sensitivity.

When pLDDT weighting or filtering is enabled, the program reads the PDB B-factor field. Interpret that field as pLDDT only when the input file actually stores prediction confidence there. Experimental temperature factors require a different interpretation.

## Pore geometry

HOLE provides radius profiles along a pore axis and the minimum radius, allowing comparison of constrictions. Keep the starting point, direction, structural inputs, and settings comparable across models.

Minimum radius is a geometric indicator. Conductance and selectivity require additional evidence. Preserve the actual HOLE inputs and logs and inspect the surrounding structure.

## Electrostatics and qualitative labels

ChimeraX scripts standardize views, surfaces, and coloring to support comparison of local electrostatics and contact changes.

`Patch_Electrostatics` and `Contacts_Qualitative` in `stage3_table.csv` summarize local surface or contact features. Check their generation or annotation rules against the code and actual outputs. Rule-based contact descriptions should be interpreted alongside counts, distances, and model quality.

## Information to retain

- Project version or Git commit.
- WT and mutant model provenance and quality information.
- Chain ID, residue numbering, residue sets, and distance cutoffs.
- Versions of ChimeraX, HOLE, and Clustal Omega used in the analysis.
- Generated scripts, logs, original metrics, and final summary tables.
- Scoring rules, cutoff-selection procedure, and experimental follow-up.

## Citation

Cite the software and model sources actually used: UCSF ChimeraX, HOLE, and the relevant structure-prediction or alignment tools. Include the repository URL and analysis commit to help others reproduce your workflow.
