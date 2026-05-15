# PCA on Breast Cancer Gene Expression Data

This project reproduces parts of Figure 1 from a breast cancer gene expression study using PCA.

The analysis uses gene expression data from 105 breast cancer patient samples to:

- Extract expression levels of **XBP1** and **GATA3**
- Visualize ER+ and ER− samples in a 2D scatter plot
- Perform PCA on the full gene-expression matrix
- Project samples onto the first principal component (PC1)

## Dataset Files

- `data/class.tsv`  
  Binary class labels:
  - `1` → ER+ breast cancer
  - `0` → ER− breast cancer

- `data/filtered.tsv.gz`  
  Gene expression matrix  
  - rows = patient samples
  - columns = gene IDs

- `data/columns.tsv.gz`  
  Mapping between gene IDs and gene symbols

## Output

The script generates:

- Scatter plot of XBP1 vs GATA3
- PCA direction visualization
- Projection of samples onto PC1

