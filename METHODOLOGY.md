# Methodology

## Overview

This document describes the detailed methodology used in each part of the scRNA-seq analysis pipeline.

---

## Part 1: Pre-processing of 10X Single-Cell RNA Datasets

### Data Source
- **Dataset**: 1k PBMCs from a Healthy Donor (v3 chemistry), 10X Genomics
- **Sub-sampled**: ~300 cells subset from original 1000 cells (Zenodo: 3457880)
- **Reference Genome**: Human hg19 chromosome X (GRCh37)
- **Whitelist**: 3M-february-2018.txt.gz (10X Genomics v3 barcodes)

### Step 1: Mapping and Demultiplexing with RNA STARsolo
- **Tool**: RNA STARsolo v2.7.11a (Galaxy)
- **Chemistry**: Chromium v3 (16bp barcode + 12bp UMI)
- **Parameters**:
  - Reference: Human hg19 chrX (built-in)
  - UMI deduplication: CellRanger2-4 algorithm
  - UMI filtering: Remove UMIs with N and homopolymers
  - Barcode matching: Multiple matches (1MM_multi)
  - Strandedness: Same strand as original RNA
  - Cell filter: Disabled (manual filtering applied later)
- **Output**: BAM alignments, count matrix (MTX), barcodes, genes, log

### Step 2: Mapping Quality Assessment with MultiQC
- **Tool**: MultiQC v1.27 (Galaxy)
- **Input**: STARsolo log file
- **Purpose**: Visualize mapping statistics and uniquely mapped read percentage

### Step 3: Cell Filtering with DropletUtils
Three filtering approaches were applied:

#### 3a. DefaultDrops (Cell Ranger method)
- Expected cells: 3000
- Upper quantile: 0.99
- Lower proportion: 0.1

#### 3b. Rank Barcodes (Knee Plot)
- Lower bound: 100
- Purpose: Identify knee and inflection points to guide threshold selection

#### 3c. EmptyDrops (Custom Filtering) — Final method used
- Lower-bound threshold: 200
- FDR threshold: 0.01
- **Result: 252 high-quality cells retained**

---

## Part 2: Basic scRNA-seq Analysis (Preprocessing & Clustering)

### Data Loading
- MTX format files loaded using `scipy.io.mmread()`
- Converted to CSR sparse matrix
- AnnData object created with barcodes as obs index and gene names as var index
- Gene file contains two columns: Ensembl ID and gene name (both retained)

### Quality Control
- **Library**: `scanpy.pp.calculate_qc_metrics()`
- **QC variables**: Ribosomal genes (RPL*, RPS* prefixes)
- **Note**: No mitochondrial genes present (chrX-only dataset)
- **Metrics computed**:
  - `n_genes_by_counts`: Number of genes detected per cell
  - `total_counts`: Total UMI counts per cell
  - `pct_counts_ribo`: Percentage of ribosomal gene counts
- **Filtering thresholds**: min_genes=50, min_cells=3 (permissive)

### Doublet Detection
- **Tool**: Scrublet (via `sc.pp.scrublet()`)
- **Method**: Nearest-neighbor classifier of observed vs simulated doublets
- **Output**: `doublet_score` and `predicted_doublet` added to obs

### Normalization
1. Raw counts saved to `adata.layers["counts"]`
2. Total count normalization to median depth: `sc.pp.normalize_total()`
3. Log1p transformation: `sc.pp.log1p()`
4. Normalized data also saved to `adata.layers["log_normalized"]`

### Feature Selection
- **Method**: Highly Variable Genes (HVG)
- **Flavor**: Seurat
- **n_top_genes**: 500 (appropriate for chrX-only dataset)

### Dimensionality Reduction
1. **PCA**: Applied on HVGs only (`use_highly_variable=True`)
2. **Neighborhood Graph**: k=15 neighbors, n_pcs=20
3. **UMAP**: Standard UMAP embedding for 2D visualization

### Clustering
- **Algorithm**: Leiden (igraph implementation, n_iterations=2)
- **Resolutions tested**: 0.1, 0.3, 0.5
- **Selected resolution**: 0.3 (best biological interpretability)

### Differential Expression & Annotation
- **Method**: Wilcoxon rank-sum test (`sc.tl.rank_genes_groups()`)
- **Groupby**: Leiden clusters
- **Annotation**: Manual based on top marker genes and known chrX immune markers

---

## Part 3: AnnData Format Exploration

### Input
- `pbmc_chrX_analyzed.h5ad` from Part 2

### Components Explored
1. **adata.X** — Active data matrix (log-normalized, sparse CSR)
2. **adata.obs** — Cell metadata DataFrame
3. **adata.var** — Gene metadata DataFrame
4. **adata.obsm** — Multi-dimensional embeddings (PCA, UMAP)
5. **adata.layers** — Raw counts, log-normalized, CPM
6. **adata.uns** — Unstructured metadata and project info
7. **adata.obsp** — Cell-cell connectivity and distance matrices

### New Layer Added: CPM Normalization
- Counts Per Million (CPM) computed from raw layer
- Formula: `CPM = (raw_count / total_cell_counts) × 1,000,000`
- Stored as `adata.layers["counts_per_million"]`

### Subsetting Demonstrated
- By numerical index
- By cell barcode name
- By boolean mask (cluster membership, QC flags)
- Combined cell + gene subsetting

### File Format
- **Format**: h5ad (HDF5-based AnnData format)
- **Compression**: gzip
- **Final output**: `pbmc_chrX_final.h5ad`

---

## Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| Galaxy | usegalaxy.org | Part 1 pipeline platform |
| RNA STARsolo | 2.7.11a | Mapping and demultiplexing |
| DropletUtils | 1.10.0 | Cell filtering |
| MultiQC | 1.27 | Quality control reporting |
| Python | 3.10+ | Parts 2 & 3 |
| Scanpy | 1.9+ | scRNA-seq analysis |
| AnnData | 0.9+ | Data structure |
| NumPy | 1.24+ | Numerical computing |
| Pandas | 2.0+ | Data manipulation |
| Scrublet | 0.2.3+ | Doublet detection |
| leidenalg | 0.10+ | Leiden clustering |
| igraph | 0.10+ | Graph algorithms |

---

## Data Flow Summary

```
[RAW INPUT]
4 FASTQ files (L001_R1, L001_R2, L002_R1, L002_R2)
Barcode whitelist (3M-february-2018.txt.gz)
Reference genome (hg19 chrX)
        ↓
[PART 1 — GALAXY]
STARsolo → 5200 detected cells
DropletUtils EmptyDrops → 252 high-quality cells
Output: barcodes.tsv + genes.tsv + matrix.mtx
        ↓
[PART 2 — SCANPY]
Load MTX → AnnData (252 cells × 2392 genes)
QC → Filter → Normalize → HVGs → PCA → UMAP → Leiden → Annotate
Output: pbmc_chrX_analyzed.h5ad
        ↓
[PART 3 — ANNDATA]
Load h5ad → Explore structure → Add CPM layer → Add project metadata
Output: pbmc_chrX_final.h5ad
```
