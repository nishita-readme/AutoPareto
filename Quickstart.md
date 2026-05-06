# AutoPareto: Quickstart

> **Copy-paste your way through for a quick tutorial.**
> For details on what each step does, see the full [README](README.md).

---

## 1. Clone & set up

```bash
git clone https://github.com/nishita-readme/AutoPareto.git
cd AutoPareto
```

Create and activate the environment:

```bash
conda env create -f environment/environment.yml
conda activate parti
```

> **No conda?**
> ```bash
> pip install scanpy partipy gseapy anndata pandas numpy matplotlib scipy
> ```

If using Notebook/Lab, register the environment as a Jupyter kernel
```bash
pip install ipykernel
python -m ipykernel install --user --name parti --display-name "Python (parti)"
```
---

## 2. Install R dependencies

In R or RStudio:

```r
install.packages(c(
  "Seurat",    # >= 5.0.0
  "Matrix",    # >= 1.6.0
  "optparse",  # >= 1.7.3
  "dplyr"      # >= 1.1.0
))
```

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "glmGamPoi",  # >= 1.14.0
  "SeuratDisk"  # >= 0.0.0.9021
))
```

---

## 3. Download test data

| Role | Accession | Files to download |
|---|---|---|
| Query | [GSE124952](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE124952) | `GSE124952_expression_matrix.csv.gz`, `GSE124952_meta_data.csv.gz` |
| Reference | [GSE115746](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE115746) | `GSE115746_exon_counts.csv.gz`, `GSE115746_complete_metadata_28706-cells.csv.gz` |

Place them in:
```
data/
├── test/raw/      ← query files
└── ref/raw/       ← reference files
```
You can do so manually, or run the following commands:
```bash
# Create directories
mkdir -p data/test/raw data/ref/raw
 
# Download query files (GSE124952)
wget -P data/test/raw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_expression_matrix.csv.gz"
wget -P data/test/raw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_meta_data.csv.gz"
 
# Download reference files (GSE115746)
wget -P data/ref/raw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746_exon_counts.csv.gz"
wget -P data/ref/raw "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115746/suppl/GSE115746_complete_metadata_28706-cells.csv.gz"
```
---

## 4. Stage 1 — Annotation Transfer (R)

> **Path note:** All file paths below assume your R working directory is `AutoPareto/notebooks/`. If it isn't, update the paths accordingly or run `setwd("path/to/AutoPareto/notebooks")` before proceeding.

> **Shortcut:** You can run all of Stage 1 at once by sourcing the tutorial script:
> ```r
> source("scripts/tutorial_stage1.R")
> ```
> Or follow the steps below to run them individually.

**Load data**

```r
source("scripts/utils_R.R")

ref <- load_seurat(
  counts_path = "data/ref/raw/GSE115746_exon_counts.csv.gz",
  meta_path   = "data/ref/raw/GSE115746_complete_metadata_28706-cells.csv.gz"
)
query <- load_seurat(
  counts_path = "data/test/raw/GSE124952_expression_matrix.csv.gz",
  meta_path   = "data/test/raw/GSE124952_meta_data.csv.gz"
)
```

**QC**

```r
ref   <- plot_qc(ref,   mt_pattern = "^mt-", ribo_pattern = "^Rp[sl]")
query <- plot_qc(query, mt_pattern = "^mt-", ribo_pattern = "^Rp[sl]")
```

```r
ref   <- subset(ref,   nFeature_RNA > 500 & pct.mt < 15)
query <- subset(query, nFeature_RNA > 500 & pct.mt < 15)
```

**Transfer annotations**

```r
query <- preprocess_and_run_transferanchor(
  query           = query,
  reference       = ref,
  normalization   = "lognorm",
  ref_label_col   = "cell_subclass",
  query_label_col = "predicted_subclass",
  dims            = 1:30,
  n_features      = 3000
)
```

**Subset to a cell type** *(optional)*

```r
query <- subset(query, subset = predicted_subclass == "L2/3 IT")
```

**Export for Python**

```r
save_annotated_data(
  seu                = query,
  counts_output_file = "data/query_counts.csv",
  meta_output_file   = "data/query_metadata.csv"
)
```

---

## 5. Stage 2 — Archetypal Analysis (Python)

> **Path note:** All file paths below assume your notebook is running from `AutoPareto/notebooks/`. If it isn't, you will need to update all `../data/...` paths to match your actual data location. You can check your current working directory with `import os; os.getcwd()`.

> **Shortcut:** You can run all of Stage 2 interactively by opening the tutorial notebook:
> ```
> notebooks/tutorial_stage2.ipynb
> ```
> Or follow the steps below to run them individually.

Before running any cells, make sure the notebook is using the `parti` kernel. In Jupyter Notebook/Lab, go to **Kernel → Change Kernel** and select **Python (parti)**. If it doesn't appear, make sure you completed the ipykernel registration step in §1.

**Load data**

```python
import scanpy as sc
import pandas as pd
import partipy as pt
import sys

# Set this to the root of your AutoPareto clone, e.g. "/home/user/AutoPareto"
REPO_ROOT = "/path/to/AutoPareto"

sys.path.append(REPO_ROOT)
from scripts.utils import *

counts = pd.read_csv("../data/query_counts.csv", index_col=0)
meta   = pd.read_csv("../data/query_metadata.csv", index_col=0)
adata  = sc.AnnData(X=counts.T, obs=meta)
```

**Verify that there are raw integer counts in adata.X**

```python
check_raw_integers_in_adataX(adata)
```

**Preprocessing**

```python
QC_genes = pd.read_csv("../data/accessories/QC_genes.txt", sep="\t").iloc[:, 0].tolist()
preprocess_adata(adata, exclude_quality_genes=True, custom_exclude_genes=QC_genes, n_pcs=50, pca_seed=123)
```

**Determine number of informative PCs**

Inspect the plot, then set `n_dims` accordingly.

```python
pt.compute_shuffled_pca(adata)
pt.plot_shuffled_pca(adata)
n_dims = 8  # adjust based on plot
```

**Select number of archetypes**

Inspect the plots, then set `n_archetypes` accordingly.

```python
pt.set_obsm(adata, obsm_key="X_pca", n_dimensions=n_dims)
pt.compute_selection_metrics(adata, n_archetypes_list=list(range(2, 8)))
plots_for_n_archetypes_selection(adata, n_archetype_range=range(3, 8), color="CellType")
n_archetypes = 3  # adjust based on plots
```

**Assign cells to archetypes**

```python
adata = get_top_cells_per_archetype(adata, n_archetypes=n_archetypes, top_n=200, n_dims=n_dims)
plot_top_cells_per_archetype(adata, dims=(0, 1))
```

**Differential expression**

```python
deg_dict          = run_deg_per_archetype(adata, lfc_threshold=1.0, pval_threshold=0.05)
pairwise_deg_dict = run_pairwise_deg_per_archetype(adata, lfc_threshold=1.0, pval_threshold=0.05)
strict_genes_df   = get_strict_archetype_genes(deg_dict, pairwise_deg_dict)
```

**GO enrichment**

```python
go_results        = run_go_analysis(deg_dict, adata, organism="mouse", n_top_genes=200)
strict_go_results = run_strict_go_analysis(strict_genes_df, adata, organism="mouse")
```

---

That's it. Results are stored in `go_results`, `strict_go_results`, `deg_dict`, and `adata`.
