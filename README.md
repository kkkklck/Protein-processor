# OsAKT2 Gating Analysis Toolkit (Python + ChimeraX)

This repository contains a Python-based graphical toolkit that automates a large fraction of the structural analysis workflow for inward-rectifier / Shaker-like K⁺ channels.  
It was originally developed for an undergraduate research project on the rice channel **OsAKT2** and its S6 gating triad mutants (e.g. DMI, DMT, GT, ND), but the workflow is general enough to be adapted to other membrane channels.

The toolkit couples a Tkinter GUI with a backend of reusable functions that:

- generate **ChimeraX `.cxc` scripts** for standardized visualization and local analysis,
- orchestrate **HOLE** calculations through WSL,
- parse **SASA** and **hydrogen-bond** logs from ChimeraX,
- merge HOLE and SASA/H-bond metrics into a single table,  
- compute simple **gating / hydration scores** and combine them with **pLDDT-based confidence**.

---

## 1. Overview

In many mutagenesis projects, most of the structural analysis is conceptually simple but technically tedious:

- prepare WT and multiple mutant PDB files,
- generate comparable views and measurements in ChimeraX,
- run HOLE on a consistent pore axis for all models,
- collect and compare SASA, hydrogen bonds and minimal radii near the gate,
- assess whether a mutant is likely to tighten or loosen the gate,  
  and whether it perturbs local hydration.

This toolkit is designed to make those steps **repeatable, scriptable and less error-prone**, while keeping every output in human-readable files (PNG, HTML, CSV) that can be inspected and reused in downstream analysis or figures.

The GUI is intentionally “small but complete”: it focuses on OsAKT2-style gating problems, yet exposes enough parameters to remain flexible for other systems.

---

## 2. Main Features

- **Research mode (WT + mutants)**
  - One-click generation of a ChimeraX `.cxc` script for:
    - loading WT and multiple mutants,
    - optional alignment,
    - Coulombic surface views,
    - local **contact** analysis around a user-defined residue set,
    - **hydrogen-bond** analysis,
    - **SASA** measurement for selected residues.
  - Automatically saves images and logs in a structured directory layout for each model.

- **Mutate mode (swapaa builder)**
  - GUI-based configuration of multiple mutations on a WT structure.
  - Generates separate **`swapaa` `.cxc` scripts** under `MUT/`, which can be executed in ChimeraX to build corresponding mutant PDBs in a reproducible way.

- **HOLE mode (WSL pipeline)**
  - Prepares a HOLE working directory from WT + mutant PDBs.
  - Writes HOLE `.inp` files with consistent geometric parameters (`cpoint`, `cvect`, `sample`, `endrad`, radius file).
  - Optionally executes HOLE through WSL and collects the resulting log files.
  - Summarizes HOLE profiles to CSV, including:
    - sample points along the pore axis,
    - minimal radius and its axial position,
    - estimated gate segment and gate length.

- **Axis helper for cpoint / cvect**
  - Builds a small ChimeraX script that:
    - focuses on selected “gate” residues,
    - computes their centre of mass as `cpoint`,
    - shifts along +Z to define `cvect`.
  - Parses the resulting ChimeraX log to auto-fill HOLE parameters, reducing manual trial-and-error.

- **SASA + H-bond summarization**
  - Parses ChimeraX logs (`*_sasa.html` and `*_hbonds.txt`),
  - Produces:
    - `sasa_hbonds_summary.csv` – per-model total SASA and H-bond counts,
    - `sasa_per_residue.csv` – per-residue SASA for the selected region,
  - Automatically reorganizes files into per-model subfolders.

- **Integrated metrics and scoring**
  - Merges HOLE and SASA/H-bond data into `metrics_all.csv`.
  - Computes:
    - **GateTightScore**: compares minimal radius and gate length to WT,
    - **HydroScore**: compares SASA and H-bonds near the gate to WT,
    - **TotalScore** = GateTightScore + HydroScore,
    - **ScoreClass**: qualitative classification (`better_than_WT`, `similar_to_WT`, `worse_than_WT`).
  - If pLDDT annotations are present and Biopython is available, summarizes:
    - mean / median pLDDT per model,
    - counts of low / medium / high-confidence residues,
    - a simple **ConfidenceClass**.
  - Writes the final table to `metrics_scored.csv`.

- **OsAKT2 MSA helper (Clustal Omega)**
  - Wraps **Clustal Omega** (native or via WSL) for a multi-sequence FASTA file.
  - Produces:
    - an alignment file (`*_OsAKT2.aln`),
    - CSV views for OsAKT2-aligned positions and automatically suggested candidate sites.

---

## 3. Code Structure

- **`PP.py`**  
  Core backend module. It implements:

  - script builders:
    - `build_cxc_script` – research mode ChimeraX script,
    - `build_mutation_cxc` – swapaa scripts for mutants,
    - `build_axis_cxc` – axis-finding script for cpoint / cvect;
  - HOLE utilities:
    - `hole_write_input`, `hole_run_in_wsl`,
    - `hole_parse_profile`, `compute_gate_metrics`,
    - `hole_summarize_logs`, `hole_plot_profiles`;
  - SASA / H-bond parsing and summarization:
    - `parse_sasa_html`, `parse_hbonds_txt`,
    - `summarize_sasa_hbonds`;
  - metrics integration:
    - `merge_all_metrics`, `score_metrics_file`,
    - optional pLDDT summarization via Biopython;
  - MSA helpers:
    - `run_osakt2_msa`, `run_osakt2_msa_wsl`.

- **`graphic(PP).py`**  
  Tkinter GUI that organises the workflow into three logical “modes”:

  - **Research** (WT + mutants, `.cxc` generation, SASA/H-bonds summarization, HOLE metrics merge and scoring),
  - **Mutate** (swapaa script generator with scrollable list of mutations),
  - **HOLE** (working directory preparation, HOLE run, plotting and axis helper).

The GUI imports and directly calls the functions in `PP.py`, so the backend can also be reused from standalone scripts or notebooks.

---

## 4. Requirements

- **Operating system**
  - Windows 10 / 11 is assumed (GUI and paths are written with Windows in mind).
  - **WSL** is used for HOLE and (optionally) Clustal Omega.

- **Python**
  - Python ≥ 3.10 (type hints use `|` unions).

- **Python packages**
  - Required:
    - `pandas`
    - `matplotlib`
  - Optional but highly recommended:
    - `biopython` (for pLDDT parsing / confidence summaries)

  Install them for example with:

  ```bash
  pip install pandas matplotlib biopython
External tools

UCSF ChimeraX
Used to execute generated .cxc scripts for visualization, SASA and H-bond logs.

HOLE
Installed in a (Conda) environment inside WSL; the toolkit calls it via bash -lc.

Clustal Omega
Either on the Windows PATH or inside WSL (e.g. /usr/bin/clustalo).

You will need to adjust the paths and environment names for HOLE and Clustal Omega in PP.py to match your local installation.

5. Installation
Clone this repository

bash
复制代码
git clone https://github.com/your-username/your-repo.git
cd your-repo
Install Python dependencies

bash
复制代码
pip install pandas matplotlib biopython
Configure HOLE and Clustal Omega

In PP.py, edit the configuration variables so that they point to:

your WSL initialization script and Conda environment (for HOLE),

the HOLE executable inside that environment,

the Clustal Omega executable (either clustalo on Windows, or /usr/bin/clustalo in WSL).

Check ChimeraX

Make sure ChimeraX is installed and that you can run:

text
复制代码
runscript path\to\your_script.cxc
inside the ChimeraX command line.

6. Usage
6.1 Launching the GUI
From the repository directory, run:

bash
复制代码
python graphic_PP.py   # replace with the actual GUI script name in this repo
A scrollable window will open with three main sections: Research, Mutate, and HOLE.

6.2 Research Mode – WT + Mutant Analysis
Inputs

WT PDB

Select the WT structure (AlphaFold or experimental) at the top of the window.

Mutant PDBs

In the mutant section, add rows with:

a label (e.g. DMI, DMT, GT, ND),

a path to the corresponding PDB file.

The labels are used consistently across output filenames and later HOLE steps.

Target residues & chain (optional but recommended)

Specify:

chain ID (e.g. A),

a residue expression such as 298,299,300 or 298-305.

These residues will be used for contacts, H-bonds and SASA analysis.

Output settings

Choose:

an output directory for images and logs,

the path to save the generated .cxc script.

Features to run

Select at least one of:

full Coulombic views,

contacts,

H-bonds,

SASA.

Generating the script

Click “Generate .cxc”.
The program writes a ChimeraX script that:

loads WT and all selected mutants,

applies the chosen analyses,

saves images and logs into the output directory with filenames derived from model labels.

Run the script in ChimeraX:

text
复制代码
runscript path\to\auto_chimerax.cxc
Wait until all measurements and images are generated.

6.3 Summarizing SASA and H-bonds
After the .cxc script has finished in ChimeraX:

Return to the GUI and keep the same output directory.

Click “Summarize SASA / H-bonds”.

The backend scans *_sasa.html and *_hbonds.txt files, then:

computes total SASA and H-bond counts per model,

extracts per-residue SASA for the selected region,

writes:

sasa_hbonds_summary.csv,

sasa_per_residue.csv,

and automatically organizes files into per-model subfolders (WT/, DMI/, DMT/, etc.), which simplifies manual inspection in ChimeraX and Excel.

6.4 HOLE Mode – Pore Geometry and Gate Metrics
Step 1 – Preparing PDBs for HOLE

Choose a HOLE working directory.

Provide the list of model names (e.g. WT,DMI,DMT,GT,ND).

Click “Prepare HOLE PDB from WT + mutants”:

WT + mutant PDBs are copied into the HOLE directory,

filenames are normalized as <Model>.pdb.

Step 2 – Axis helper (optional but recommended)

Specify chain and gate residue expression (e.g. A, 298-300).

Click “Generate axis script”:

a small .cxc script is written to the HOLE directory.

Run that script in ChimeraX; it will log the estimated cpoint and cvect to a log file.

Use the parsed values as cpoint/cvect in the HOLE parameter fields.

Step 3 – Running HOLE

Set:

cpoint and cvect,

sampling step (sample),

ending radius (endrad),

radius file (e.g. simple.rad),

HOLE command (inside WSL).

Tick the checkboxes to:

write HOLE input files,

optionally execute HOLE,

parse logs.

For each model, a subdirectory <Model>-HOLE is created with:

the copied PDB,

HOLE input file(s),

HOLE log(s).

HOLE profiles are parsed and summarized into CSV tables:

hole_profile_samples.csv

hole_min_table.csv

hole_min_summary.csv

The GUI can additionally generate basic comparison plots (e.g. overlaid radius profiles and simple bar plots for minimal radius / gate length) directly in the HOLE directory.

6.5 Merging HOLE and SASA Metrics and Scoring
Once HOLE results are available and SASA/H-bond summary files have been generated:

Ensure that:

HOLE directory is set (for HOLE CSVs),

output directory in Research mode is set (for SASA/H-bond CSVs).

Click “Merge HOLE + SASA metrics”.

This will:

read hole_min_table.csv and sasa_hbonds_summary.csv,

merge them into metrics_all.csv under the output directory,

run score_metrics_file to compute:

GateTightScore,

HydroScore,

TotalScore,

ScoreClass,

optionally augment the table with pLDDT-based statistics (if PDB files with pLDDT fields are available and Biopython is installed),

finally write metrics_scored.csv.

metrics_scored.csv is intended to be the main table for ranking and selecting mutants, and can be imported directly into Excel, R or Python notebooks.

6.6 Mutate Mode – Swapaa Script Builder
In the Mutate tab:

Set a WT PDB and an output directory.

Add one or more mutation rows, each with:

a label (e.g. DMI),

one or more chain IDs,

one or more residue indices,

one or more target amino acids.

The GUI supports comma- or space-separated lists (with matching lengths) to define multiple positions in a single row.

Click the button to generate swapaa scripts:

.cxc files are written into <output_dir>/MUT/,

a summary dialog lists the generated paths and reminds you that in ChimeraX you can simply run:

text
复制代码
runscript path\to\some_mutation.cxc
Each script builds the corresponding mutant model from the WT structure, saving a PDB file with a name derived from the label.

6.7 OsAKT2 MSA Helper
For OsAKT2-related sequence analyses, the MSA helper uses Clustal Omega to:

align a multi-FASTA file of AKT/AKT-like sequences,

generate a .aln file (*_OsAKT2.aln),

call an external OsAKT2-specific tool to output:

an alignment view CSV (alignment_osakt2_view.csv),

an automatically curated list of candidate positions (candidate_sites_auto_v0.1.csv).

This module is specific to the underlying OsAKT2 project but can be adapted if you have similar alignment post-processing scripts.

7. Typical Workflow (OsAKT2 S6 Gating Mutants)
A typical use case for this toolkit looks like:

Build mutant PDBs

Use Mutate mode to generate swapaa scripts and construct DMI, DMT, GT, ND (etc.) from the WT structure in ChimeraX.

Run standardized structural analysis

Use Research mode to generate a .cxc script for WT + mutants.

Run the script in ChimeraX to obtain Coulombic views, contacts, H-bonds and SASA logs.

Summarize SASA and H-bonds

Use “Summarize SASA / H-bonds” to generate CSV summaries and per-residue SASA tables.

Analyze pore geometry with HOLE

Use HOLE mode to prepare PDBs, get cpoint / cvect, run HOLE through WSL, and generate radius profiles.

Merge metrics and score mutants

Use “Merge HOLE + SASA metrics” to obtain metrics_all.csv and metrics_scored.csv.

Inspect and rank mutants by TotalScore, comparing gating tightness, hydration and structural confidence to WT.

Document and visualize

Use the generated PNGs and CSVs as the basis for figures and tables in reports or manuscripts.

8. Acknowledgements
This toolkit is a thin automation layer around several powerful open-source tools and resources:

UCSF ChimeraX – visualization, SASA measurements and hydrogen-bond analysis.

HOLE – pore radius profiling.

Clustal Omega – multiple sequence alignment.

AlphaFold / PDB – structural inputs (where applicable).

When using this toolkit in scientific work, please cite the underlying tools according to their respective documentation and publications.
