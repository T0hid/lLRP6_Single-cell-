## 📦 Data Availability & Input Files

### Raw Data Sources
This analysis uses scRNA-seq data from tooth development stages. The raw count matrices can be downloaded from NCBI GEO:

| Stage | Sample ID | Description | Link |
|-------|-----------|-------------|------|
| E13.5 | `GSM5700357` | Lower Molar Mesenchyme | [View on NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5700357) |
| E14.5 | `GSM5700358` | Lower Molar Mesenchyme | [View on NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5700358) |
| E16.5 | `GSM5700359` | Lower Molar Mesenchyme | [View on NCBI](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5700359) |

### Processed Files
To run the final visualization script (`lig-rep.py`), ensure the following local files are present:

* **Processed Object:** `processed_data/processed_adata.h5ad`
* **Interaction Scores:** [`interactions_FINAL/lrp5_lrp6_interactions.csv`](./interactions_FINAL/lrp5_lrp6_interactions.csv)
