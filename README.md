# Protein Pipeline (PP) — Reproducible Structure-Analysis Workflow for Point Mutants

**Protein Pipeline (PP)** is a GUI-driven, reproducible workflow for **batch structural analysis of point mutations**, originally developed for ion channels (e.g., Shaker-type K⁺ channels) but designed to be **protein-agnostic** as long as you provide PDB models with consistent residue numbering.

It helps you go from “WT + mutant structures” → to **standardized figures + quantitative metrics + summary tables**, so downstream wet-lab validation has clear, testable hypotheses.

---

## What PP Can Do (Current Capabilities)

### 1) Generate standardized UCSF ChimeraX scripts (`.cxc`)
PP generates ChimeraX command scripts to produce **consistent, comparable outputs** across WT and mutants, such as:

- **Local electrostatics around a site/patch** (e.g., P-loop pocket): surface Coulombic coloring and site-focused views
- **Local contact visualization**: residue sticks + local contact networks around user-defined residues
- **SASA / H-bond summaries** for defined regions (site-centric or residue-list centric)
- **Export images** (`sites_coulombic_img`, `sites_contacts_img`, etc.) in a consistent style for reporting

> The goal is not only visualization, but **standardization** (same viewpoint + same settings) so differences are interpretable.

---

### 2) Pore geometry analysis with HOLE (via WSL)
PP runs **HOLE** in a Windows + WSL setup and summarizes results into:
- `hole_min_table.csv` (per-model minimum radius, position, etc.)
- `hole_min_summary.csv` (compact summary for quick comparison)
- optional plots (hole profile curves, bar charts)

This provides a **direct geometric proxy** for conductance-relevant constraints (e.g., narrowing near the gate/selectivity region).

---

### 3) “Cross contacts” quantification (P-loop ↔ S6)
PP includes a contact quantification module intended to replace vague “contacts look dense/loose” language with **numbers**.

Typical outputs:
- `CrossContactPairs`: count of residue pairs across two residue groups passing a distance cutoff
- `CrossContactDensity`: normalized count (e.g., pairs / theoretical maximum pairs)
- `CrossContactMinDist`: minimum inter-group distance observed
- optional cutoff sweep columns: `Pairs@4.0`, `Pairs@5.0`, `Pairs@6.0`, `Pairs@6.5`, `Pairs@7.0`

Exports:
- `contacts_cross_summary.csv`
- optionally merged into `metrics_all.csv` / `stage3_table.csv`

---

### 4) Project-level summary tables (for reporting & scoring)
PP aggregates outputs into:
- `metrics_all.csv` (all computed metrics in one place)
- `metrics_scored.csv` (metrics mapped into a scoring system if enabled)
- `stage3_table.csv` (human-readable reporting table, including qualitative tags)

---

## What Those “Stage 3” Qualitative Columns Mean (No Code Needed)

In `stage3_table.csv`, two columns translate raw structural outputs into plain-language evidence.

### `Patch_Electrostatics`
A qualitative description of the **electrostatic patch** around a defined pocket/region (e.g., P-loop neighborhood), typically read from the **Coulombic surface images**.

Example tags:
- “moderately negative + half-open”
- “strongly negative + concave/half-open”

Interpretation:
> Indicates whether the pocket tends to attract/repel charged species and whether its shape appears open/occluded.

---

### `Contacts_Qualitative`
A qualitative description of whether the **coupling contact network** between two functional regions is intact, guided by contact images and supported by cross-contact metrics.

Example tags:
- “contacts weakened”
- “contacts broken”
- “contacts compact/intact”

Interpretation:
> Summarizes whether the P-loop ↔ S6 coupling appears mechanically connected or decoupled.

---

## Methods-Style Note: Cross Contacts Metrics and Cutoff Sweep (How to Interpret Results)

### Definitions
Given two residue sets **A** (e.g., P-loop residues) and **B** (e.g., S6 residues), PP computes inter-residue proximity using atomic coordinates from PDB models.

- For each residue pair *(i ∈ A, j ∈ B)*, PP computes a distance metric (commonly the **minimum heavy-atom distance** between the two residues, or a chosen representative atom distance depending on the implementation).
- A residue pair is counted as a “cross contact” if its distance is ≤ **cutoff** (Å).

From these pairwise checks, PP reports:
- **CrossContactPairs**: number of contacting residue pairs across the two sets
- **CrossContactDensity**: `CrossContactPairs / (|A| × |B|)` (normalization to allow comparison across different residue-set sizes)
- **CrossContactMinDist**: the minimum distance observed across all cross-set pairs (useful as a “closest approach” indicator)

### Why values may look identical across models
Cross contacts can appear identical for multiple models for practical reasons:
1. **Residue sets are small**: if |A| and |B| are small, `|A|×|B|` is small and **counts become coarse** (e.g., 0/1/2…), limiting resolution.
2. **Cutoff is not discriminative**:  
   - Too strict → most models report 0 pairs  
   - Too loose → most models report the same saturated pair count  
3. **Distance metric choice**: using Cα–Cα distances can miss side-chain-specific differences; heavy-atom or side-chain distances are usually more sensitive.

### Cutoff sweep as a robustness check (recommended)
To avoid overfitting a single cutoff, PP can compute a **cutoff sweep**:
- `Pairs@4.0`, `Pairs@5.0`, `Pairs@6.0`, `Pairs@6.5`, `Pairs@7.0`, …

Interpretation:
- If models diverge only at larger cutoffs (e.g., ≥6.5 Å), that suggests **weak/long-range proximity** rather than tight packing.
- If divergence exists at smaller cutoffs (e.g., 4–5 Å), that suggests **strong, physically tight coupling**.

Recommended reporting practice:
- Use `CrossContactMinDist` + a **small set of cutoff points** (e.g., 5.0, 6.5, 7.0 Å) to show where separation emerges.
- Treat cross-contact metrics as **supporting evidence** that complements gate geometry (HOLE) and electrostatics, rather than as a standalone “proof”.

---

## Repository Structure (Typical)

- `graphic（PP）.py`  
  Main GUI entry point (Tkinter). Launch this for daily use.

- `PP.py`  
  Core backend: generating ChimeraX scripts, running HOLE via WSL, extracting/summarizing metrics, table export.

- `help_texts.py`  
  In-app help/manual text.

- `msa_consensus_tool.py`  
  Utility for MSA consensus / mutation suggestion support (optional).

---

## Requirements

### Mandatory
- **Python 3.9+** (Windows recommended)
- Typical Python packages:
  - `pandas`, `numpy`
  - `matplotlib` (if plotting is enabled)
  - Tkinter (usually bundled with standard Python on Windows)

### External tools (module-dependent)
- **UCSF ChimeraX** (required for `.cxc` execution and figure generation)
- **WSL2** (Windows Subsystem for Linux)
- **HOLE** installed inside WSL (required for pore analysis)
- Optional: **Clustal Omega** (for MSA utilities)

---

## Installation (Typical)

1. Create a Python environment (conda or venv).
2. Install dependencies:
   ```bash
   pip install pandas numpy matplotlib
Install and verify UCSF ChimeraX.

If you need HOLE:

enable WSL2

install HOLE in WSL

ensure PP’s WSL paths/env activation match your setup

Quick Start (Recommended Workflow)
Step 1 — Prepare input models
Put WT + mutants in one folder, with:

consistent residue numbering

consistent chain IDs (or known chain ID to use)

Step 2 — Launch PP GUI
bash
复制代码
python "graphic（PP）.py"
Step 3 — Generate ChimeraX scripts and images
Configure in GUI:

PDB folder

chain ID

site residues (e.g., P-loop residues)

optional paired region residues (e.g., P-loop group vs S6 group)
Generate .cxc scripts and run them in ChimeraX to produce:

sites_coulombic_img

sites_contacts_img

SASA/H-bond outputs (if enabled)

Step 4 — Run HOLE (optional)
Batch-run pore analysis and export:

hole_min_table.csv

hole_min_summary.csv

plots (if enabled)

Step 5 — Export metrics tables
Export/merge:

metrics_all.csv

metrics_scored.csv (if enabled)

stage3_table.csv

Reproducibility & Batch Scaling
PP is designed so a structural workflow can be replayed and scaled:

standardized .cxc scripts instead of manual clicking

batch outputs per model

merged CSV summaries for ranking and reporting

persistent UI state (e.g., settings.json) to avoid reconfiguring each session

This makes it feasible to go from “6 mutants” → “60 mutants” without turning the project into chaos.

Troubleshooting (Common Failure Modes)
Residue not found / numbering mismatch:
Most issues come from inconsistent residue numbering across models.

Chain mismatch:
Different chain IDs across models can break site extraction and contact analysis.

HOLE/WSL failures:
Usually path/env activation issues inside WSL.

Cross-contact metrics all identical:
Indicates a non-discriminative cutoff, a too-small residue set, or an insensitive distance metric.
Use cutoff sweep columns and report CrossContactMinDist to show where separation emerges.

Scientific Note (How to Use Outputs Responsibly)
PP does not “prove” functional change alone. It provides:

standardized structural evidence,

quantitative ranking signals,

mechanistic narratives (electrostatics + geometry + coupling contacts),

so wet-lab assays can validate the most plausible candidates first.

How to Cite (Guideline)
If you use PP in academic work, cite the underlying tools:

UCSF ChimeraX (structure visualization/analysis)

HOLE (pore radius profiling)

your structure predictor used for input models (e.g., AlphaFold / ColabFold)

License
This project is released under the MIT License.

Contact / Contribution
The pipeline evolves alongside real research workflows. Issues and PRs are welcome, especially for:

more robust residue mapping across models

improved contact metrics (heavy-atom-only, sidechain-only, per-residue min-dist, etc.)

data-driven cutoff selection and better sensitivity analysis
