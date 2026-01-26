# scRNA-seq Analysis Pipeline with scVI Integration

A complete single-cell RNA sequencing analysis workflow built around scVI for batch correction and imputation. Takes raw 10X Genomics data through QC, clustering, cell type annotation, and ligand-receptor interaction analysis.

## What it does

This pipeline handles the full journey from raw count matrices to publication-ready figures:

- Quality control filtering (gene counts, UMI counts, mitochondrial content)
- Batch correction and denoising using scVI deep learning
- Dimensionality reduction and UMAP visualisation
- Leiden clustering with automated cell type annotation based on marker genes
- Ligand-receptor interaction scoring
- Statistical validation using both imputed and raw count data

The code is set up for Wnt signalling analysis but can be adapted to any pathway or gene set of interest.

## Requirements

Python 3.8+ with:

```bash
pip install scanpy scvi-tools matplotlib seaborn pandas numpy scipy
```

GPU recommended for scVI training but not required.

## Data setup

Organise your 10X Genomics outputs by sample/condition:

```
processed_stages/
├── Sample1/
│   ├── matrix.mtx.gz
│   ├── features.tsv.gz
│   └── barcodes.tsv.gz
├── Sample2/
│   └── ...
└── Sample3/
    └── ...
```

Update the `SAMPLE_MAP` dictionary in the `Config` class to match your sample names.

## Running

```bash
python scrnaseq_pipeline.py
```

Takes roughly 30-45 minutes depending on dataset size and hardware.

## Output

Results go to `final_manuscript_output/`:

- `figures/` — QC plots, marker validation, UMAPs, and expression summaries
- `data/final_processed_revised.h5ad` — complete AnnData object with all annotations and embeddings

## Customisation

Edit the `Config` class to adjust:

- **QC thresholds**: `MIN_GENES`, `MIN_CELLS`, `MAX_MT`
- **scVI parameters**: layer count, latent dimensions, training epochs
- **Clustering**: Leiden resolution, neighbour count
- **Gene lists**: update `WNT_LIGANDS`, `WNT_RECEPTORS` for your pathway
- **Cell type markers**: modify the `MARKERS` dictionary for your tissue

## Statistical approach

The pipeline computes expression statistics on multiple data representations: scVI-imputed values, standard log-normalised counts, and depth-downsampled counts. This helps distinguish genuine biological signal from technical artefacts like variable sequencing depth across cells.

## License

MIT

## Contact

Open an issue for questions or bugs.
