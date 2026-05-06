source("../scripts/utils_R.R")

ref <- load_seurat(
  counts_path = "data/ref/raw/GSE115746_exon_counts.csv.gz",
  meta_path   = "data/ref/raw/GSE115746_complete_metadata_28706-cells.csv.gz",
)

query <- load_seurat(
  counts_path = "data/test/raw/GSE124952_expression_matrix.csv.gz",
  meta_path   = "data/test/raw/GSE124952_meta_data.csv.gz",
)

## QC: adjust as needed
ref   <- plot_qc(ref,   mt_pattern = "^mt-", ribo_pattern = "^Rp[sl]")
query <- plot_qc(query, mt_pattern = "^mt-", ribo_pattern = "^Rp[sl]")

ref   <- subset(ref,   nFeature_RNA > 500 & pct.mt < 15)
query <- subset(query, nFeature_RNA > 500 & pct.mt < 15)

## Annotation Transfer
query <- preprocess_and_run_transferanchor(
  query           = query,
  reference       = ref,
  normalization   = "lognorm",        # or "SCT"
  ref_label_col   = "cell_subclass",  # column in ref metadata to transfer
  query_label_col = "predicted_subclass",
  dims            = 1:30,
  n_features      = 3000
)

## For this tutorial, subset for L2/3 IT
query <- subset(query, subset = predicted_subclass == "L2/3 IT")

## Save subsetted and annotated counts and metadata matrices
save_annotated_data(
  seu                = query,
  counts_output_file = "data/query_counts.csv",
  meta_output_file   = "data/query_metadata.csv"
)
