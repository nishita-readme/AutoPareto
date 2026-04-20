source("../scripts/utils_R.R")

ref <- load_seurat(
  counts_path = "data/ref/raw/counts.csv.gz",
  meta_path   = "data/ref/raw/metadata.csv.gz",
)

query <- load_seurat(
  counts_path = "data/test/query_counts.csv.gz",
  meta_path   = "data/test/query_metadata.csv.gz",
)

query <- preprocess_and_run_transferanchor(
  query           = query,
  reference       = ref,
  normalization   = "lognorm",        # or "SCT"
  ref_label_col   = "cell_subclass",  # column in ref metadata to transfer
  query_label_col = "predicted_subclass",
  dims            = 1:30,
  n_features      = 3000
)

save_annotated_data(
  seu                = query,
  counts_output_file = "data/processed/query_counts.csv",
  meta_output_file   = "data/processed/query_metadata.csv"
)
