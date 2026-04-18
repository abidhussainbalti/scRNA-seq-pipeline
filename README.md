# scRNA-seq Analysis Pipeline: From Raw Data to AnnData

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://www.python.org/)
[![Scanpy](https://img.shields.io/badge/Scanpy-1.9+-green.svg)](https://scanpy.readthedocs.io/)
[![Galaxy](https://img.shields.io/badge/Galaxy-usegalaxy.org-orange.svg)](https://usegalaxy.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

## Overview

This repository presents a **complete end-to-end single-cell RNA sequencing (scRNA-seq) analysis pipeline** using **1,000 Peripheral Blood Mononuclear Cells (PBMCs)** from a healthy human donor, sequenced using the **10X Genomics Chromium v3** platform.

The pipeline is divided into three interconnected parts — the output of each part feeds directly into the next, forming a coherent and reproducible analysis workflow:

```
Part 1: Raw FASTQ Files (10X Genomics)
        ↓  [RNA STARsolo + DropletUtils on Galaxy]
        ↓  Output: barcodes.tsv + genes.tsv + matrix.mtx
        ↓
Part 2: Count Matrix → Quality Control → Clustering → Cell Type Annotation
        ↓  [Scanpy in Python]
        ↓  Output: pbmc_chrX_analyzed.h5ad
        ↓
Part 3: AnnData Format Exploration → Enriched Metadata → Final Dataset
           [AnnData in Python]
           Output: pbmc_chrX_final.h5ad
```

---

## Repository Structure

```
scRNA-seq-project/
│
├── README.md                          ← You are here
├── METHODOLOGY.md                     ← Detailed methods description
├── requirements.txt                   ← Python dependencies
├── .gitignore                         ← Git ignore rules
│
├── 01-10X-Preprocessing/              ← PART 1: Galaxy Pipeline
│   ├── README.md                      ← Part 1 documentation
│   ├── workflow/
│   │   └── scrna_preprocessing.ga     ← Exported Galaxy workflow
│   └── results/
│       ├── barcodes.tsv               ← Cell barcodes (252 cells)
│       ├── genes.tsv                  ← Gene IDs and names (2392 genes)
│       └── matrix.mtx                 ← Sparse count matrix
│
├── 02-Scanpy-Clustering/              ← PART 2: Scanpy Analysis
│   ├── README.md                      ← Part 2 documentation
│   ├── part2_scanpy_clustering.ipynb  ← Main analysis notebook
│   └── results/
│       ├── pbmc_chrX_analyzed.h5ad    ← Fully analyzed AnnData object
│       └── pbmc_cell_annotations.csv  ← Cell metadata as CSV
│
├── 03-AnnData-Tutorial/               ← PART 3: AnnData Exploration
│   ├── README.md                      ← Part 3 documentation
│   ├── part3_anndata_tutorial.ipynb   ← AnnData tutorial notebook
│   └── results/
│       └── pbmc_chrX_final.h5ad       ← Final enriched AnnData object
│
└── assets/
    └── pipeline_overview.png          ← Pipeline diagram
```

---

## Part 1: Pre-processing of 10X Single-Cell RNA Datasets

**Platform**: [Galaxy](https://usegalaxy.org/) | **Tools**: RNA STARsolo, DropletUtils, MultiQC

### What was done
Raw FASTQ files from 10X Genomics were processed to produce a clean, filtered count matrix ready for downstream analysis.

### Dataset
| Property | Value |
|----------|-------|
| Sample | 1k PBMCs from Healthy Donor |
| Chemistry | 10X Chromium v3 |
| Genome | Human hg19 (chrX) |
| Sequencing lanes | L001, L002 |
| Input files | 4 FASTQ files (R1 + R2 per lane) |

### Tools & Steps
| Step | Tool | Purpose |
|------|------|---------|
| Mapping & Demultiplexing | RNA STARsolo | Align reads to genome, assign to cells |
| Quality Control | MultiQC | Mapping quality assessment |
| Cell Filtering | DropletUtils (DefaultDrops) | Cell Ranger-style filtering |
| Barcode Ranking | DropletUtils (Rank Barcodes) | Knee plot for cell quality |
| Custom Filtering | DropletUtils (EmptyDrops) | Custom threshold filtering |

### Key Results
| Metric | Value |
|--------|-------|
| Total cells detected (STARsolo) | 5,200 |
| High-quality cells (EmptyDrops) | **252** |
| Genes quantified | **2,392** (chrX) |
| Non-zero matrix entries | 39,002 |
| Filtering method | EmptyDrops (lower=200, FDR=0.01) |

### Output Files
- `results/barcodes.tsv` — 252 cell barcodes
- `results/genes.tsv` — 2,392 gene IDs and names
- `results/matrix.mtx` — Sparse count matrix (MatrixMarket format)

📁 [Go to Part 1 →](01-10X-Preprocessing/)

---

## Part 2: Basic scRNA-seq Analysis (Preprocessing & Clustering)

**Platform**: Google Colab / Jupyter | **Language**: Python | **Library**: Scanpy

### What was done
The count matrix from Part 1 was loaded into an AnnData object and processed through a complete single-cell analysis workflow.

### Pipeline Steps
| Step | Method | Details |
|------|--------|---------|
| Data Loading | scipy.io + AnnData | MTX format → AnnData |
| Quality Control | `sc.pp.calculate_qc_metrics` | Ribosomal gene % analysis |
| Cell/Gene Filtering | `sc.pp.filter_cells/genes` | min_genes=50, min_cells=3 |
| Doublet Detection | Scrublet | Predict and flag doublets |
| Normalization | `sc.pp.normalize_total` | Median count normalization |
| Log Transform | `sc.pp.log1p` | Variance stabilization |
| Feature Selection | `sc.pp.highly_variable_genes` | Top 500 HVGs |
| Dimensionality Reduction | PCA + UMAP | 50 PCs, 15 neighbors |
| Clustering | Leiden algorithm | Resolutions: 0.1, 0.3, 0.5 |
| Marker Genes | Wilcoxon rank-sum test | Per-cluster DEGs |
| Cell Type Annotation | Manual (marker-guided) | Broad PBMC cell types |

### Key Results
| Metric | Value |
|--------|-------|
| Input cells | 252 |
| Final cells after filtering | 252 |
| Highly variable genes | 500 |
| PCA components used | 20 |
| Number of clusters (res=0.3) | varies |
| Doublets detected | < 5% |

### Output Files
- `part2_scanpy_clustering.ipynb` — Complete analysis notebook with outputs
- `results/pbmc_chrX_analyzed.h5ad` — Analyzed AnnData (→ input for Part 3)
- `results/pbmc_cell_annotations.csv` — Cell metadata table

📁 [Go to Part 2 →](02-Scanpy-Clustering/)

---

## Part 3: Getting Started with AnnData

**Platform**: Google Colab / Jupyter | **Language**: Python | **Library**: AnnData

### What was done
The analyzed AnnData object from Part 2 was used to explore and demonstrate all components of the AnnData format — the standard data structure for single-cell analysis in the Python ecosystem.

### Components Explored
| Component | Description | Our Data |
|-----------|-------------|----------|
| `adata.X` | Active data matrix | Log-normalized counts (252 × 2392) |
| `adata.obs` | Cell annotations | QC metrics, clusters, cell types |
| `adata.var` | Gene annotations | Ensembl IDs, HVG flags |
| `adata.obsm` | Cell embeddings | X_pca (50D), X_umap (2D) |
| `adata.layers` | Data versions | Raw counts, log-norm, CPM |
| `adata.uns` | Unstructured metadata | Project info, analysis params |
| `adata.obsp` | Cell-cell matrices | k-NN connectivities, distances |

### Key Concepts Demonstrated
- Loading and inspecting h5ad files
- Subsetting by index, name, and boolean masks
- Views vs copies — memory efficiency
- Adding new layers (CPM normalization)
- Exporting to DataFrames and CSV
- Writing enriched AnnData back to disk

### Output Files
- `part3_anndata_tutorial.ipynb` — Complete tutorial notebook with outputs
- `results/pbmc_chrX_final.h5ad` — Final enriched AnnData object

📁 [Go to Part 3 →](03-AnnData-Tutorial/)

---

## How to Reproduce

### Requirements
```bash
pip install -r requirements.txt
```

### Run Part 2 (Scanpy Clustering)
1. Open `02-Scanpy-Clustering/part2_scanpy_clustering.ipynb` in Jupyter or Google Colab
2. Upload the 3 files from `01-10X-Preprocessing/results/` when prompted
3. Run all cells
4. Download `pbmc_chrX_analyzed.h5ad`

### Run Part 3 (AnnData Tutorial)
1. Open `03-AnnData-Tutorial/part3_anndata_tutorial.ipynb` in Jupyter or Google Colab
2. Upload `pbmc_chrX_analyzed.h5ad` from Part 2 when prompted
3. Run all cells

---

## References

1. **10X Genomics** — PBMC dataset: [1k PBMCs from a Healthy Donor (v3 chemistry)](https://www.10xgenomics.com/resources/datasets)
2. **Galaxy Training Network** — [Pre-processing of 10X Single-Cell RNA Datasets](https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-preprocessing-tenx/tutorial.html)
3. **scverse** — [Preprocessing and Clustering Tutorial](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/clustering.html)
4. **AnnData** — [Getting Started with AnnData](https://anndata.readthedocs.io/en/latest/tutorials/notebooks/getting-started.html)
5. **scverse** — [AnnData Getting Started](https://scverse-tutorials.readthedocs.io/en/latest/notebooks/anndata_getting_started.html)
6. Dobin A. et al. (2013) STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*
7. Lun A. et al. (2019) EmptyDrops: distinguishing cells from empty droplets. *Genome Biology*
8. Wolf F.A. et al. (2018) SCANPY: large-scale single-cell gene expression data analysis. *Genome Biology*

---

## License
This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

## Author
Abid Hussain-Special Topics in Bioinformatics 

