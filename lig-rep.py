#!/usr/bin/env python3
"""

Author: Tohid Ghasemnejad
Date: 2026
"""

import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import numpy as np
import os
import warnings
import scvi
from scipy import stats
from scipy.sparse import issparse

# =============================================================================
# CONFIGURATION
# =============================================================================

class Config:
    # --- PATHS ---
    BASE_PATH = os.getcwd() 
    # Ensure these folders exist or point to your specific data location
    INPUT_PATH = os.path.join(BASE_PATH, "processed_stages")
    OUTPUT_PATH = os.path.join(BASE_PATH, "final_manuscript_output")
    
    # --- QC PARAMS ---
    MIN_GENES = 200
    MIN_CELLS = 3
    MAX_MT = 20
    
    # --- scVI PARAMS ---
    SCVI_LAYERS = 2
    SCVI_LATENT = 30
    SCVI_EPOCHS = 400
    BATCH_KEY = "Stage"
    
    # --- ANALYSIS PARAMS ---
    LEIDEN_RES = 0.5
    N_NEIGHBORS = 30
    
    # --- GENE LISTS ---
    WNT_RECEPTORS = ['Lrp5', 'Lrp6']
    WNT_LIGANDS = ['Wnt3', 'Wnt4', 'Wnt5a', 'Wnt6', 'Wnt10a', 'Wnt11']
    
    # Marker Dictionary for Annotation
    MARKERS = {
        'Epithelial': ['Krt14', 'Epcam', 'Cdh1'],
        'Mesenchyme': ['Col1a1', 'Vim', 'Pdgfra', 'Twist1'],
        'Endothelial': ['Pecam1', 'Cdh5', 'Kdr'],
        'Fibroblasts': ['Col3a1', 'Dcn', 'Lum'],
        'Chondrogenic': ['Sox9', 'Col2a1', 'Acan'],
        'Osteogenic': ['Runx2', 'Sp7', 'Alpl'],
        'Macrophages': ['Cd68', 'Adgre1', 'Csf1r'],
        'Glial': ['Sox10', 'Mbp', 'Plp1']
    }
    
    # Stage Mapping
    SAMPLE_MAP = {
        'GSM5700357': 'E13.5',
        'GSM5700358': 'E14.5', 
        'GSM5700359': 'E16.5'
    }

    @staticmethod
    def setup_dirs():
        """Creates output directories if they don't exist."""
        os.makedirs(Config.OUTPUT_PATH, exist_ok=True)
        os.makedirs(os.path.join(Config.OUTPUT_PATH, "figures"), exist_ok=True)
        os.makedirs(os.path.join(Config.OUTPUT_PATH, "data"), exist_ok=True)
        
    @staticmethod
    def setup_plotting():
        """Configures matplotlib for publication-quality figures."""
        warnings.filterwarnings('ignore')
        sc.settings.verbosity = 0
        sc.settings.set_figure_params(dpi=300, facecolor='white', frameon=False)
        plt.rcParams['font.family'] = 'DejaVu Sans'
        plt.rcParams['svg.fonttype'] = 'none' # Editable text in Illustrator

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def safe_expression(adata_obj, gene, layer=None):
    """
    Safely extract flattened expression array for a gene.
    Handles sparse matrices, dense arrays, and missing genes.
    """
    if gene not in adata_obj.var_names:
        return np.array([])
    
    if layer:
        X = adata_obj[:, gene].layers[layer]
    else:
        X = adata_obj[:, gene].X
        
    if issparse(X):
        X = X.todense()
    return np.asarray(X).flatten()

# =============================================================================
# PIPELINE STEPS
# =============================================================================

def step1_load_qc(config):
    print("--- Step 1: Loading & QC ---")
    adatas = []
    
    # Load individual samples
    for gsm, stage in config.SAMPLE_MAP.items():
        try:
            p = os.path.join(config.INPUT_PATH, stage)
            # cache=True speeds up subsequent re-runs
            a = sc.read_10x_mtx(p, var_names='gene_symbols', cache=True)
            a.var_names_make_unique()
            a.obs['Stage'] = stage
            adatas.append(a)
            print(f"Loaded {stage}: {a.n_obs} cells")
        except Exception as e:
            print(f"Warning: Skipping {stage} ({e})")
    
    if not adatas:
        raise ValueError("No data found. Check your INPUT_PATH.")
    
    # Concatenate
    adata = sc.concat(adatas, label="batch_indices")
    adata.var_names_make_unique()
    
    # Calculate QC metrics
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    
    # Plot Figure S1 (QC) BEFORE filtering
    print("Generating Figure S1 (QC)...")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for i, metric in enumerate(['n_genes_by_counts', 'total_counts', 'pct_counts_mt']):
        sc.pl.violin(adata, metric, groupby='Stage', ax=axes[i], show=False, palette='viridis')
    plt.tight_layout()
    plt.savefig(os.path.join(config.OUTPUT_PATH, "figures", "FigureS1_QC.png"))
    plt.close()

    # Filter
    print(f"Pre-filter: {adata.n_obs} cells")
    sc.pp.filter_cells(adata, min_genes=config.MIN_GENES)
    adata = adata[adata.obs.pct_counts_mt < config.MAX_MT].copy()
    sc.pp.filter_genes(adata, min_cells=config.MIN_CELLS)
    
    # Save raw counts for scVI and rigorous stats
    adata.layers["counts"] = adata.X.copy() 
    print(f"Post-QC: {adata.n_obs} cells")
    return adata

def step2_scvi(adata, config):
    print("--- Step 2: scVI Deep Learning Imputation ---")
    
    # Setup AnnData for scVI
    # We use a copy for setup to avoid modifying the main object in unexpected ways
    # but scvi modifies in place mostly.
    
    # Log-normalize for standard visualization (saved to .X)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata # Freeze raw normalized state
    
    # Select Highly Variable Genes (HVGs) based on counts
    sc.pp.highly_variable_genes(
        adata, 
        n_top_genes=3000, 
        layer="counts", 
        flavor="seurat_v3", 
        batch_key=config.BATCH_KEY
    )
    
    # Train Model
    scvi.model.SCVI.setup_anndata(
        adata, 
        layer="counts", 
        batch_key=config.BATCH_KEY, 
        continuous_covariate_keys=["pct_counts_mt"]
    )
    
    model = scvi.model.SCVI(adata, n_layers=config.SCVI_LAYERS, n_latent=config.SCVI_LATENT)
    
    print("Training scVI model (this may take time)...")
    # Use GPU if available, else CPU
    model.train(
        max_epochs=config.SCVI_EPOCHS, 
        accelerator="gpu", 
        devices=1, 
        early_stopping=True
    )
    
    # Extract Latent Representation (for UMAP/Clustering)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    
    # Extract Imputed Expression (denoised)
    # We save this to a layer so we don't overwrite the raw data
    adata.layers["imputed"] = model.get_normalized_expression(library_size=1e4)
    
    # Compute Neighbors & UMAP on scVI latent space
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=config.N_NEIGHBORS)
    sc.tl.umap(adata)
    
    return adata

def step3_cluster_annotate(adata, config):
    print("--- Step 3: Clustering & Annotation ---")
    sc.tl.leiden(adata, resolution=config.LEIDEN_RES)
    
    # Automated Annotation based on Marker Score
    cluster_map = {}
    print("Annotating clusters...")
    
    for cl in adata.obs['leiden'].unique():
        subset = adata[adata.obs['leiden'] == cl]
        
        # Calculate mean expression of markers for each cell type
        scores = {}
        for ct, markers in config.MARKERS.items():
            valid_markers = [m for m in markers if m in adata.var_names]
            if not valid_markers:
                scores[ct] = 0
                continue
            # Use raw normalized data for annotation to be safe
            mean_expr = np.mean([safe_expression(subset, m).mean() for m in valid_markers])
            scores[ct] = mean_expr
            
        best_ct = max(scores, key=scores.get)
        # Threshold to avoid forcing noise into a category
        cluster_map[cl] = best_ct if scores[best_ct] > 0.1 else "Unknown"
    
    adata.obs['Cell_Type'] = adata.obs['leiden'].map(cluster_map)
    
    # Remove Unknowns for cleaner plots
    adata = adata[adata.obs['Cell_Type'] != "Unknown"].copy()
    
    # Generate Figure S2 (Markers Validation)
    print("Generating Figure S2 (Markers)...")
    flat_markers = [m for sublist in config.MARKERS.values() for m in sublist]
    # Filter markers present in data
    flat_markers = [m for m in flat_markers if m in adata.var_names]
    
    sc.pl.dotplot(adata, flat_markers, groupby='Cell_Type', standard_scale='var', show=False)
    plt.savefig(os.path.join(config.OUTPUT_PATH, "figures", "FigureS2_Markers.png"), bbox_inches='tight')
    plt.close()
    
    return adata

def step4_interaction_stats(adata, config):
    print("--- Step 4: Interaction Analysis ---")
    interactions = []
    cell_types = adata.obs['Cell_Type'].unique()
    
    for lig in config.WNT_LIGANDS:
        if lig not in adata.var_names: continue
        
        for rec in config.WNT_RECEPTORS:
            for source in cell_types:
                for target in cell_types:
                    s_data = adata[adata.obs['Cell_Type'] == source]
                    t_data = adata[adata.obs['Cell_Type'] == target]
                    
                    # Use IMPUTED data here for cleaner means (standard practice for interaction potential)
                    l_expr = safe_expression(s_data, lig, "imputed").mean()
                    r_expr = safe_expression(t_data, rec, "imputed").mean()
                    
                    if l_expr > 0.1 and r_expr > 0.1: # Threshold
                        # Geometric mean of Ligand(Source) * Receptor(Target)
                        score = np.sqrt(l_expr * r_expr)
                        interactions.append({
                            'Ligand': lig, 'Receptor': rec, 
                            'Source': source, 'Target': target, 'Score': score
                        })
    return pd.DataFrame(interactions)

def step5_rigorous_stats(adata, config):
    print("--- Step 5: Rigorous Statistical Controls (Reviewer Response) ---")
    
    # Focus on Mesenchyme (the relevant lineage)
    mes_cells = adata[adata.obs['Cell_Type'] == 'Mesenchyme'].copy()
    print(f"Analyzing {mes_cells.n_obs} Mesenchymal cells...")
    
    # --- 1. RAW DATA ANALYSIS ---
    # We must use raw counts, re-normalize, and log-transform 
    # WITHOUT using scVI smoothing to get the "honest" correlation.
    mes_raw = mes_cells.copy()
    mes_raw.X = mes_raw.layers['counts'].copy() 
    sc.pp.normalize_total(mes_raw, target_sum=1e4)
    sc.pp.log1p(mes_raw)
    
    x_raw = safe_expression(mes_raw, 'Lrp6')
    y_raw = safe_expression(mes_raw, 'Lrp5')
    r_raw, p_raw = stats.pearsonr(x_raw, y_raw)
    
    # --- 2. DOWNSAMPLING CONTROL ---
    # Address concern: "Is correlation just driven by read depth?"
    # We downsample all cells to the same read depth (25th percentile)
    target_counts = int(np.percentile(mes_cells.obs['total_counts'], 25))
    
    mes_ds = mes_cells.copy()
    mes_ds.X = mes_ds.layers['counts'].copy()
    sc.pp.downsample_counts(mes_ds, total_counts=target_counts, replace=False)
    sc.pp.log1p(mes_ds) # Log transform after downsampling
    
    x_ds = safe_expression(mes_ds, 'Lrp6')
    y_ds = safe_expression(mes_ds, 'Lrp5')
    r_ds, p_ds = stats.pearsonr(x_ds, y_ds)

    print(f"\n>>> STATISTICAL SUMMARY:")
    print(f"1. Raw Data (Log-Norm): r = {r_raw:.3f} (p={p_raw:.2e})")
    print(f"2. Downsampled Control: r = {r_ds:.3f} (p={p_ds:.2e})")
    
    return {
        'x_raw': x_raw, 'y_raw': y_raw, 'r_raw': r_raw, 'p_raw': p_raw,
        'x_ds': x_ds, 'y_ds': y_ds, 'r_ds': r_ds, 'p_ds': p_ds
    }

# =============================================================================
# FIGURE GENERATION
# =============================================================================

def generate_figure_4(adata, config):
    print("Generating Figure 4 (Spatiotemporal)...")
    fig = plt.figure(figsize=(18, 15))
    gs = gridspec.GridSpec(3, 2, width_ratios=[1, 1.2])
    
    stages = ['E13.5', 'E14.5', 'E16.5']
    genes = ['Lrp6', 'Lrp5', 'Wnt10a', 'Axin2', 'Msx1']
    
    for i, stage in enumerate(stages):
        subset = adata[adata.obs['Stage'] == stage]
        
        # UMAP Panel
        ax_umap = fig.add_subplot(gs[i, 0])
        if subset.n_obs > 0:
            sc.pl.umap(subset, color='Cell_Type', ax=ax_umap, show=False, legend_loc='on data', 
                      frameon=False, title=f"{stage} - Cell Types")
        
        # Dotplot Panel
        ax_dot = fig.add_subplot(gs[i, 1])
        if subset.n_obs > 0:
            sc.pl.dotplot(subset, genes, groupby='Cell_Type', ax=ax_dot, show=False, 
                         title=f"{stage} - Wnt Components", cmap='Reds', standard_scale='var')
                      
    plt.tight_layout()
    plt.savefig(os.path.join(config.OUTPUT_PATH, "figures", "Figure4_Spatiotemporal.png"))
    plt.close()

def generate_figure_5_revised(adata, interactions_df, stats_data, config):
    """
    Generates the composite figure for the manuscript.
    Includes the 'Reviewer Control' panel (Raw Data) side-by-side with Imputed.
    """
    print("Generating Revised Figure 5 (Response to Reviewers)...")
    fig = plt.figure(figsize=(24, 12)) 
    gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1]) 
    
    # --- PANEL A1: Imputed Co-expression (The "Model" View) ---
    ax_imp = fig.add_subplot(gs[0, 0])
    mes = adata[adata.obs['Cell_Type'] == 'Mesenchyme']
    x_imp = safe_expression(mes, 'Lrp6', 'imputed')
    y_imp = safe_expression(mes, 'Lrp5', 'imputed')
    
    # Calculate stats for imputed just for the label
    r_imp, p_imp = stats.pearsonr(x_imp, y_imp)
    
    hb1 = ax_imp.hexbin(x_imp, y_imp, gridsize=40, cmap='Blues', mincnt=1, bins='log')
    m, b = np.polyfit(x_imp, y_imp, 1)
    ax_imp.plot(x_imp, m*x_imp + b, 'r--', lw=2)
    ax_imp.set_title(f"A1. Imputed Expression (scVI)\nModel-Denoised\nr={r_imp:.2f}")
    ax_imp.set_xlabel("Lrp6 (Imputed)")
    ax_imp.set_ylabel("Lrp5 (Imputed)")
    plt.colorbar(hb1, ax=ax_imp, label='Log10(Count)')

    # --- PANEL A2: Raw Co-expression (The "Rigorous" View) ---
    # 
    ax_raw = fig.add_subplot(gs[0, 1])
    x_raw, y_raw = stats_data['x_raw'], stats_data['y_raw']
    r_raw, p_raw = stats_data['r_raw'], stats_data['p_raw']
    
    hb2 = ax_raw.hexbin(x_raw, y_raw, gridsize=40, cmap='Greens', mincnt=1, bins='log')
    if len(x_raw) > 0:
        m2, b2 = np.polyfit(x_raw, y_raw, 1)
        ax_raw.plot(x_raw, m2*x_raw + b2, 'k--', lw=2)
        
    ax_raw.set_title(f"A2. Raw Normalized Counts\nReviewer Control\nr={r_raw:.2f}, p={p_raw:.2e}")
    ax_raw.set_xlabel("Lrp6 (Log-Normalized)")
    ax_raw.set_ylabel("Lrp5 (Log-Normalized)")
    plt.colorbar(hb2, ax=ax_raw, label='Log10(Count)')

    # --- PANEL B: Receptor Dynamics ---
    ax_dyn = fig.add_subplot(gs[0, 2])
    data = []
    for s in ['E13.5', 'E14.5', 'E16.5']:
        sub = adata[adata.obs['Stage'] == s]
        # Use RAW data for bar plots to be conservative/accurate to depth
        sub_raw = sub.copy()
        sub_raw.X = sub_raw.layers['counts'].copy()
        sc.pp.normalize_total(sub_raw, target_sum=1e4)
        sc.pp.log1p(sub_raw)
        
        data.append({'Stage': s, 'Gene': 'Lrp5', 'Exp': safe_expression(sub_raw, 'Lrp5').mean()})
        data.append({'Stage': s, 'Gene': 'Lrp6', 'Exp': safe_expression(sub_raw, 'Lrp6').mean()})
    
    if data:
        sns.barplot(data=pd.DataFrame(data), x='Stage', y='Exp', hue='Gene', ax=ax_dyn, palette=['#3498db', '#e74c3c'])
        ax_dyn.set_title("B. Receptor Dynamics (Raw)")
        ax_dyn.set_ylabel("Mean Log-Norm Expression")

    # --- PANEL C: Ligand Heatmap ---
    ax_heat = fig.add_subplot(gs[0, 3])
    lig_mat = pd.DataFrame(index=adata.obs['Cell_Type'].unique(), columns=config.WNT_LIGANDS)
    for ct in lig_mat.index:
        sub = adata[adata.obs['Cell_Type'] == ct]
        for lig in config.WNT_LIGANDS:
            if lig in adata.var_names:
                # Use imputed means for heatmap clarity (standard practice)
                lig_mat.loc[ct, lig] = safe_expression(sub, lig, 'imputed').mean()
    
    sns.heatmap(lig_mat.astype(float), cmap="Reds", ax=ax_heat, cbar=False)
    ax_heat.set_title("C. Ligand Source")

    # --- PANELS D & E: Interactions ---
    ax_int6 = fig.add_subplot(gs[1, :2])
    if not interactions_df.empty:
        df6 = interactions_df[interactions_df['Receptor'] == 'Lrp6'].sort_values('Score', ascending=False).head(10)
        df6['Label'] = df6['Ligand'] + " (" + df6['Source'].str[:4] + "->" + df6['Target'].str[:4] + ")"
        sns.barplot(data=df6, y='Label', x='Score', ax=ax_int6, color='#e74c3c')
        ax_int6.set_title("D. Top LRP6 Interactions")

    ax_int5 = fig.add_subplot(gs[1, 2:])
    if not interactions_df.empty:
        df5 = interactions_df[interactions_df['Receptor'] == 'Lrp5'].sort_values('Score', ascending=False).head(10)
        df5['Label'] = df5['Ligand'] + " (" + df5['Source'].str[:4] + "->" + df5['Target'].str[:4] + ")"
        sns.barplot(data=df5, y='Label', x='Score', ax=ax_int5, color='#3498db')
        ax_int5.set_title("E. Top LRP5 Interactions")

    plt.tight_layout()
    plt.savefig(os.path.join(config.OUTPUT_PATH, "figures", "Figure5_Revised_Mechanism.png"))
    plt.close()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    Config.setup_dirs()
    Config.setup_plotting()
    print(">>> STARTING PIPELINE")
    
    # 1. Load & QC
    adata = step1_load_qc(Config)
    
    # 2. scVI Imputation (Generates 'imputed' layer and 'X_scVI' latent)
    adata = step2_scvi(adata, Config)
    
    # 3. Cluster & Annotate
    adata = step3_cluster_annotate(adata, Config)
    
    # 4. Interactions
    int_df = step4_interaction_stats(adata, Config)
    
    # 5. NEW: Run Rigorous Statistics
    stats_results = step5_rigorous_stats(adata, Config)
    
    # 6. Generate Figures
    generate_figure_4(adata, Config)
    
    # Pass the stats results to the figure generator for Panel A2
    generate_figure_5_revised(adata, int_df, stats_results, Config)
    
    # Save Data
    print("Saving processed object...")
    adata.write(os.path.join(Config.OUTPUT_PATH, "data", "final_processed_revised.h5ad"))
    
    print(">>> PIPELINE COMPLETE.")
    print("Check 'final_manuscript_output/figures' for revised plots.")
