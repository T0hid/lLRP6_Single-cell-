#!/usr/bin/env python3
"""
IMPROVED FIGURE: Shows LRP5 and LRP6 interactions SEPARATELY
=============================================================

"""

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

# Try to import adjustText, install if not available
try:
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True
except ImportError:
    HAS_ADJUSTTEXT = False
    print("Note: adjustText not installed. Installing...")
    import subprocess
    subprocess.check_call(['pip', 'install', 'adjustText', '-q'])
    from adjustText import adjust_text
    HAS_ADJUSTTEXT = True

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10


def create_improved_figure(adata, lrp_interactions, output_dir):
    """
    Create figure with separate LRP5 and LRP6 interaction panels.
    Panel F (compensation model) has been removed.
    """
    
    # Adjusted figure size for 2 rows instead of 3
    fig = plt.figure(figsize=(18, 11))
    gs = fig.add_gridspec(2, 3, hspace=0.35, wspace=0.4,
                          height_ratios=[1, 1.2])
    
    cell_types = sorted(adata.obs['Cell_Type'].unique())
    
    # =========================================================================
    # PANEL A: Co-expression (with adjustText for non-overlapping labels)
    # =========================================================================
    ax1 = fig.add_subplot(gs[0, 0])
    
    expr_data = []
    for ct in cell_types:
        subset = adata[adata.obs['Cell_Type'] == ct]
        if 'Lrp6' in adata.var_names and 'Lrp5' in adata.var_names:
            expr_data.append({
                'Cell_Type': ct,
                'Lrp6': float(subset[:, 'Lrp6'].X.mean()),
                'Lrp5': float(subset[:, 'Lrp5'].X.mean())
            })
    
    df_expr = pd.DataFrame(expr_data)
    
    # Plot scatter points
    ax1.scatter(df_expr['Lrp6'], df_expr['Lrp5'], 
               s=120, c='steelblue', alpha=0.7, edgecolors='white', linewidth=1.5,
               zorder=10)
    
    # Expand axis limits BEFORE adding labels to give room
    x_data = df_expr['Lrp6'].values
    y_data = df_expr['Lrp5'].values
    
    x_min, x_max = 0, max(x_data) * 1.25
    y_min, y_max = 0, max(y_data) * 1.15
    ax1.set_xlim(x_min, x_max)
    ax1.set_ylim(y_min, y_max)
    
    # Create text objects for adjustText
    texts = []
    for _, row in df_expr.iterrows():
        txt = ax1.text(row['Lrp6'], row['Lrp5'], row['Cell_Type'],
                      fontsize=8, alpha=0.9, zorder=5)
        texts.append(txt)
    
    # Use adjustText to prevent overlaps
    adjust_text(texts, 
                x=df_expr['Lrp6'].values,
                y=df_expr['Lrp5'].values,
                ax=ax1,
                arrowprops=dict(arrowstyle='-', color='gray', lw=0.5, alpha=0.5),
                expand_points=(1.5, 1.5),
                expand_text=(1.2, 1.2),
                force_points=(0.5, 0.5),
                force_text=(0.8, 0.8),
                only_move={'points': 'y', 'texts': 'xy'},
                ensure_inside_axes=True,
                time_lim=3)
    
    # Add correlation coefficient
    corr = df_expr['Lrp6'].corr(df_expr['Lrp5'])
    ax1.text(0.05, 0.95, f'r = {corr:.2f}', transform=ax1.transAxes,
            fontsize=11, fontweight='bold',
            bbox=dict(facecolor='white', edgecolor='gray', alpha=0.9))
    
    ax1.set_xlabel('Lrp6 Expression')
    ax1.set_ylabel('Lrp5 Expression')
    ax1.set_title('A. LRP5/LRP6 Co-expression', fontweight='bold')
    
    # =========================================================================
    # PANEL B: Receptor dynamics (same as before)
    # =========================================================================
    ax2 = fig.add_subplot(gs[0, 1])
    
    stage_expr = []
    for stage in ['E13.5', 'E14.5', 'E16.5']:
        subset = adata[adata.obs['Stage'] == stage]
        for gene in ['Lrp5', 'Lrp6']:
            if gene in adata.var_names:
                stage_expr.append({
                    'Stage': stage, 
                    'Gene': gene, 
                    'Expression': float(subset[:, gene].X.mean())
                })
    
    df_stage = pd.DataFrame(stage_expr)
    x = np.arange(3)
    width = 0.35
    
    for i, (gene, color) in enumerate([('Lrp5', '#3498db'), ('Lrp6', '#e74c3c')]):
        data = df_stage[df_stage['Gene'] == gene]['Expression'].values
        if len(data) > 0:
            ax2.bar(x + i*width, data, width, label=gene, color=color, alpha=0.85,
                   edgecolor='white', linewidth=0.5)
    
    ax2.set_xlabel('Developmental Stage')
    ax2.set_ylabel('Mean Expression')
    ax2.set_title('B. Receptor Expression Dynamics', fontweight='bold')
    ax2.set_xticks(x + width/2)
    ax2.set_xticklabels(['E13.5', 'E14.5', 'E16.5'])
    ax2.legend(frameon=True, loc='upper right')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    
    # =========================================================================
    # PANEL C: Wnt ligand expression (same as before)
    # =========================================================================
    ax3 = fig.add_subplot(gs[0, 2])
    
    wnt_genes = ['Wnt10a', 'Wnt11', 'Wnt3', 'Wnt4', 'Wnt5a', 'Wnt6']
    available_wnt = [g for g in wnt_genes if g in adata.var_names]
    
    if available_wnt:
        wnt_expr = []
        for ct in cell_types:
            subset = adata[adata.obs['Cell_Type'] == ct]
            for gene in available_wnt:
                wnt_expr.append({
                    'Cell_Type': ct,
                    'Gene': gene,
                    'Expression': float(subset[:, gene].X.mean())
                })
        
        df_wnt = pd.DataFrame(wnt_expr)
        df_pivot = df_wnt.pivot(index='Cell_Type', columns='Gene', values='Expression')
        df_pivot = df_pivot[available_wnt]
        
        im = ax3.imshow(df_pivot.values, cmap='Reds', aspect='auto')
        ax3.set_xticks(np.arange(len(available_wnt)))
        ax3.set_yticks(np.arange(len(cell_types)))
        ax3.set_xticklabels(available_wnt, rotation=45, ha='right', fontsize=9)
        ax3.set_yticklabels(cell_types, fontsize=9)
        
        # Improved colorbar positioning
        cbar = plt.colorbar(im, ax=ax3, shrink=0.8, pad=0.02)
        cbar.set_label('Expression', fontsize=9)
    
    ax3.set_title('C. Wnt Ligand Expression', fontweight='bold')
    
    # =========================================================================
    # PANEL D: SIDE-BY-SIDE LRP6 vs LRP5 INTERACTIONS
    # =========================================================================
    
    if lrp_interactions is not None and len(lrp_interactions) > 0:
        # Split by receptor
        lrp6_int = lrp_interactions[lrp_interactions['receptor'] == 'Lrp6'].head(10)
        lrp5_int = lrp_interactions[lrp_interactions['receptor'] == 'Lrp5'].head(10)
        
        # Panel D: LRP6 interactions
        ax4a = fig.add_subplot(gs[1, 0:2])
        
        if len(lrp6_int) > 0:
            y_pos = np.arange(len(lrp6_int))
            scores = lrp6_int['interaction_score'].values / lrp_interactions['interaction_score'].max()
            
            bars = ax4a.barh(y_pos, scores, color='#e74c3c', alpha=0.85, 
                            edgecolor='white', linewidth=0.5)
            labels = [f"{row['ligand']} ({row['source']} → {row['target']})"
                     for _, row in lrp6_int.iterrows()]
            ax4a.set_yticks(y_pos)
            ax4a.set_yticklabels(labels, fontsize=9)
            ax4a.set_xlabel('Interaction Score (normalized)')
            ax4a.set_title('D. LRP6 Interactions (Primary Receptor)', fontweight='bold', color='#e74c3c')
            ax4a.set_xlim(0, 1.1)
            ax4a.spines['top'].set_visible(False)
            ax4a.spines['right'].set_visible(False)
            
            # Add value labels on bars
            for bar, score in zip(bars, scores):
                if score > 0.1:
                    ax4a.text(score + 0.02, bar.get_y() + bar.get_height()/2,
                             f'{score:.2f}', va='center', fontsize=8, color='gray')
        
        # Panel E: LRP5 interactions
        ax4b = fig.add_subplot(gs[1, 2])
        
        if len(lrp5_int) > 0:
            y_pos = np.arange(len(lrp5_int))
            # Normalize to same scale as LRP6 for fair comparison
            scores = lrp5_int['interaction_score'].values / lrp_interactions['interaction_score'].max()
            
            bars = ax4b.barh(y_pos, scores, color='#3498db', alpha=0.85, 
                            edgecolor='white', linewidth=0.5)
            labels = [f"{row['ligand']} ({row['source'][:4]}→{row['target'][:4]})"
                     for _, row in lrp5_int.iterrows()]
            ax4b.set_yticks(y_pos)
            ax4b.set_yticklabels(labels, fontsize=8)
            ax4b.set_xlabel('Interaction Score')
            ax4b.set_title('E. LRP5 Interactions\n(Compensatory)', fontweight='bold', color='#3498db')
            ax4b.set_xlim(0, 1.1)
            ax4b.spines['top'].set_visible(False)
            ax4b.spines['right'].set_visible(False)
            
            # Add note about lower scores
            ax4b.text(0.5, -0.12, 'Note: Lower scores reflect\nlower Lrp5 expression',
                     transform=ax4b.transAxes, fontsize=8, ha='center',
                     style='italic', color='gray')
            
            # Add value labels on bars
            for bar, score in zip(bars, scores):
                if score > 0.05:
                    ax4b.text(score + 0.02, bar.get_y() + bar.get_height()/2,
                             f'{score:.2f}', va='center', fontsize=7, color='gray')
        
        # Find shared ligands
        lrp6_ligands = set(lrp6_int['ligand'].unique())
        lrp5_ligands = set(lrp5_int['ligand'].unique())
        shared = lrp6_ligands.intersection(lrp5_ligands)
        
        if shared:
            # Add annotation about shared ligands
            ax4a.text(0.02, 0.98, f'Shared with LRP5: {", ".join(sorted(shared))}',
                     transform=ax4a.transAxes, fontsize=9, va='top',
                     bbox=dict(facecolor='#ffffcc', edgecolor='orange', alpha=0.8))
    else:
        # If no interaction data, show placeholder
        ax4a = fig.add_subplot(gs[1, 0:2])
        ax4a.text(0.5, 0.5, 'No interaction data available\nRun interaction analysis first',
                 ha='center', va='center', transform=ax4a.transAxes,
                 fontsize=12, color='gray')
        ax4a.axis('off')
        
        ax4b = fig.add_subplot(gs[1, 2])
        ax4b.axis('off')
    
    # Main title
    fig.suptitle('Single-cell Evidence for LRP5/LRP6 Compensation in Tooth Development',
                fontsize=14, fontweight='bold', y=0.98)
    
    # Save
    os.makedirs(output_dir, exist_ok=True)
    
    png_path = os.path.join(output_dir, 'Figure_LRP_Compensation_v2.png')
    plt.savefig(png_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {png_path}")
    
    pdf_path = png_path.replace('.png', '.pdf')
    plt.savefig(pdf_path, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"Saved: {pdf_path}")
    
    plt.close()
    
    # Print summary for manuscript
    if lrp_interactions is not None:
        lrp6_ligands = set(lrp_interactions[lrp_interactions['receptor'] == 'Lrp6']['ligand'].head(10))
        lrp5_ligands = set(lrp_interactions[lrp_interactions['receptor'] == 'Lrp5']['ligand'].head(10))
        shared = sorted(lrp6_ligands.intersection(lrp5_ligands))
        
        print("\n" + "="*60)
        print("FOR YOUR ABSTRACT/RESULTS:")
        print("="*60)
        print(f"\nShared ligands targeting BOTH LRP6 and LRP5: {shared}")
        print(f"\nSuggested text:")
        print(f'"Ligand-receptor interaction analysis revealed that {", ".join(shared)}')
        print('can signal through both LRP6 and LRP5 in odontogenic cell populations,')
        print('providing a mechanistic basis for LRP5-mediated compensation."')


# =============================================================================
# MAIN
# =============================================================================

if __name__ == "__main__":
    base_path = os.getcwd()
    
    # Load data
    h5ad_path = os.path.join(base_path, 'processed_data', 'processed_adata.h5ad')
    int_path = os.path.join(base_path, 'interactions_FINAL', 'lrp5_lrp6_interactions.csv')
    
    if not os.path.exists(h5ad_path):
        print(f"ERROR: Data not found at {h5ad_path}")
        print("Run lrp6_analysis_WORKING.py first!")
        exit(1)
    
    print("Loading data...")
    adata = sc.read_h5ad(h5ad_path)
    
    if os.path.exists(int_path):
        lrp_interactions = pd.read_csv(int_path)
        print(f"Loaded {len(lrp_interactions)} interactions")
    else:
        # Calculate interactions if not available
        print("Calculating interactions...")
        # ... (you'd need to include the calculation code here)
        lrp_interactions = None
    
    output_dir = os.path.join(base_path, 'figures_FINAL')
    create_improved_figure(adata, lrp_interactions, output_dir)
    
    print("\nDone!")