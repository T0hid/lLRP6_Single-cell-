
Single-Cell Wnt Signaling & LRP5/LRP6 Interaction Analysis
This repository contains Python pipelines for analyzing single-cell RNA sequencing (scRNA-seq) data in tooth development. The code focuses on visualizing gene expression trends across developmental stages (E13.5, E14.5, E16.5) and quantifying LRP5/LRP6 ligand-receptor interactions to investigate compensatory mechanisms.

📂 Repository Contents
1. multi_gene_panel.py (formerly better-vis.py)
Function: End-to-end preprocessing and multi-stage visualization.

Workflow:

Ingests raw 10x count matrices for E13.5, E14.5, and E16.5.

Performs QC, normalization, log-transformation, and clustering (Leiden resolution 0.4).

Renames clusters to biological annotations (e.g., Mesenchyme, Epithelial, Odontoblasts) and standardizes colors across stages.

Output: Generates a composite figure with UMAPs and Dot Plots for target genes (Lrp6, Lrp5, Wnt10a, Axin2, Msx1).

2. lig-rep.py (formerly lrp6-final.py)
Function: Generates the final publication-quality figure for Ligand-Receptor analysis.

Workflow:

Panel A: Scatter plot of Lrp6/Lrp5 co-expression with non-overlapping labels (via adjustText).

Panel B: Bar charts showing receptor expression dynamics across stages.

Panel C: Heatmap of Wnt ligand expression across cell types.

Panel D/E: Comparative bar charts of interaction scores, highlighting shared ligands between LRP6 (primary) and LRP5 (compensatory).

Output: Saves a high-resolution figure and prints a text summary of shared ligands for manuscript use.

🛠️ Installation
The scripts require Python 3 and the following dependencies:

Bash

pip install scanpy matplotlib pandas numpy adjustText
Note: lig-rep.py attempts to auto-install adjustText if it is missing, but manual installation is recommended.

📁 Directory Structure
Ensure your directory is structured as follows for the scripts to locate input files correctly:

Plaintext

.
├── multi_gene_panel.py          <-- Preprocessing & Basic Vis
├── lig-rep.py                   <-- Final Interaction Figure
├── processed_stages/            <-- INPUT for multi_gene_panel.py
│   ├── E13.5/
│   ├── E14.5/
│   └── E16.5/
├── processed_data/              <-- INPUT for lig-rep.py
│   └── processed_adata.h5ad     <-- Pre-calculated AnnData object
└── interactions_FINAL/          <-- INPUT for lig-rep.py
    └── lrp5_lrp6_interactions.csv
🚀 Usage
Step 1: Preprocessing & Gene Panels
Run multi_gene_panel.py to process the raw stage data and generate the initial visualization panels.

Bash

python multi_gene_panel.py
Output: figures_LRP_MultiGene/final_multi_gene_panel.png

Customization: Edit target_genes in the script to visualize different markers.

Step 2: Final Interaction Figure
Run lig-rep.py to generate the detailed LRP5 vs LRP6 comparison figure.

Bash

python lig-rep.py
Output: figures_FINAL/Figure_LRP_Compensation_v2.png (and PDF).

Console Output: The script will print a summary of shared ligands (e.g., "Shared ligands targeting BOTH LRP6 and LRP5") to assist in writing results.

⚙️ Configuration Notes
Cluster Renaming: multi_gene_panel.py contains a dictionary cluster_renaming mapping numeric clusters to cell types. Modify this if your clustering resolution changes.

Interaction Data: lig-rep.py expects a CSV containing columns for receptor, ligand, source, target, and interaction_score. If this file is missing, the script will output a placeholder figure for the interaction panels.
