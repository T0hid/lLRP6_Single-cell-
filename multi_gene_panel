import scanpy as sc
import matplotlib.pyplot as plt
import os

# --- CONFIGURATION ---
target_genes = ['Lrp6', 'Lrp5', 'Wnt10a', 'Axin2', 'Msx1']
# 1. FIX: Added '15' and '10' (if needed) to ensure clean filtering
clusters_to_remove = ['11', '12', '13', '14', '15', 'Unknown']

base_path = os.getcwd()
input_path = os.path.join(base_path, "processed_stages")
figures_folder = os.path.join(base_path, "figures_LRP_MultiGene")
os.makedirs(figures_folder, exist_ok=True)

# --- LOAD AND PREPROCESS ---
# (Your loading code remains the same)
adatas = []
sample_map = {'GSM5700357': 'E13.5', 'GSM5700358': 'E14.5', 'GSM5700359': 'E16.5'}

for gsm_id, stage in sample_map.items():
    stage_dir = os.path.join(input_path, stage)
    try:
        adata = sc.read_10x_mtx(stage_dir, var_names='gene_symbols', cache=True)
        adata.var_names_make_unique()
        adata.obs['Stage'] = stage
        adatas.append(adata)
    except Exception as e:
        print(f"Skipping {stage}: {e}")

adata_atlas = sc.concat(adatas, label="batch", keys=sample_map.values())
adata_atlas.obs_names_make_unique()

# Standard processing
sc.pp.filter_cells(adata_atlas, min_genes=200)
sc.pp.filter_genes(adata_atlas, min_cells=3)
sc.pp.normalize_total(adata_atlas, target_sum=1e4)
sc.pp.log1p(adata_atlas)
sc.pp.highly_variable_genes(adata_atlas, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.pca(adata_atlas)
sc.pp.neighbors(adata_atlas)
sc.tl.umap(adata_atlas)
sc.tl.leiden(adata_atlas, resolution=0.4) 

# --- RENAME CLUSTERS ---
# 2. FIX: Ensure all numeric clusters are accounted for
cluster_renaming = {
    '0': 'Mesenchyme', '1': 'Epithelial', '2': 'Fibroblasts',
    '3': 'Endothelial', '4': 'Myogenic', '5': 'Chondrogenic',
    '6': 'Osteogenic', '7': 'Neural Crest', '8': 'Macrophages',
    '9': 'Unknown', '10': 'Glial', '11': 'Unknown', '12': 'Unknown',
    '13': 'Unknown', '14': 'Unknown', '15': 'Unknown' 
}

# Apply renaming
current_clusters = adata_atlas.obs['leiden'].unique()
safe_renaming = {k: v for k, v in cluster_renaming.items() if k in current_clusters}
adata_atlas.obs['Cell_Type'] = adata_atlas.obs['leiden'].map(safe_renaming).fillna(adata_atlas.obs['leiden'])

# Filter
adata_atlas = adata_atlas[~adata_atlas.obs['Cell_Type'].isin(clusters_to_remove)].copy()

# --- 3. FIX: CONSISTENT COLORS ---
# Create a fixed color map so E13.5 Mesenchyme looks the same as E16.5 Mesenchyme
unique_types = sorted(adata_atlas.obs['Cell_Type'].unique())
# Use the 'tab20' palette for distinct colors
palette = plt.get_cmap('tab20')
color_map = {ctype: palette(i) for i, ctype in enumerate(unique_types)}
adata_atlas.uns['Cell_Type_colors'] = [color_map[c] for c in unique_types] # Pre-load for scanpy

# --- PLOTTING (CORRECTED) ---
stages = ['E13.5', 'E14.5', 'E16.5']

# 1. Increase figure width slightly and adjust spacing
fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(18, 14), gridspec_kw={'width_ratios': [1, 0.7]})
plt.subplots_adjust(wspace=0.25, hspace=0.4) # Add breathing room

valid_genes = [g for g in target_genes if g in adata_atlas.var_names]

for i, stage in enumerate(stages):
    subset = adata_atlas[adata_atlas.obs['Stage'] == stage].copy()
    
    # Ensure consistent colors for this subset
    subset.uns['Cell_Type_colors'] = [color_map[c] for c in sorted(subset.obs['Cell_Type'].unique())]
    
    # --- LEFT PANEL: UMAP ---
    # CHANGE: legend_loc='on data' puts labels directly on the clusters
    sc.pl.umap(
        subset, 
        color='Cell_Type', 
        ax=axes[i, 0], 
        show=False, 
        title=f"{stage} - Clusters",
        frameon=False, 
        legend_loc='on data',      # <--- Key Fix: Places text on the clusters
        legend_fontsize=9,         # slightly smaller text
        legend_fontoutline=2,      # white outline to make text readable against color
        palette=color_map, 
        size=20                    # dots slightly larger
    )
    
    # --- RIGHT PANEL: DOT PLOT ---
    # The dot plot will now have clear space to the left
    sc.pl.dotplot(
        subset, 
        valid_genes, 
        groupby='Cell_Type', 
        ax=axes[i, 1],
        show=False,
        title=f"{stage} - Expression",
        vmin=0, vmax=1.5,
        cmap='Reds',
        colorbar_title='Mean Exp.'
    )
    
    # Clean up axis labels
    axes[i, 1].set_ylabel('') 

output_file = os.path.join(figures_folder, "final_multi_gene_panel.png")
plt.savefig(output_file, dpi=300, bbox_inches='tight')
print(f"\nSUCCESS: Figure saved to {output_file}")
