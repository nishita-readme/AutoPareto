library(Seurat)
library(data.table)
library(ggplot2)
library(patchwork)

#' @title Load a counts matrix and optional metadata into a Seurat object
#'
#' @description
#' Reads a counts matrix and optional metadata file from disk and creates
#' a Seurat object. Cell IDs are matched between counts columns and metadata
#' rownames automatically.
#'
#' @param counts_path Path to counts/expression matrix (csv, csv.gz, tsv, tsv.gz).
#'   Expected orientation: rows = genes, cols = cells.
#' @param meta_path Optional path to metadata file. First column is used as cell
#'   IDs and must overlap with column names of the counts matrix.
#' @param transpose Logical. Set TRUE if matrix is cells x genes instead of
#'   genes x cells. Default FALSE.
#' @param project Project name passed to CreateSeuratObject, stored in orig.ident.
#'   Default "SeuratProject".
#'
#' @return A Seurat object.
#'
#' @examples
#' seu <- load_seurat(
#'   counts_path = "data/counts.csv.gz",
#'   meta_path   = "data/metadata.csv.gz",
#'   project     = "MyExperiment"
#' )
load_seurat <- function(counts_path,
                        meta_path = NULL,
                        transpose = FALSE,
                        project   = "SeuratProject") { ... }


#' @title Plot QC metrics for a Seurat object
#'
#' @description
#' Computes mitochondrial and ribosomal percentages per cell and produces two
#' rows of diagnostic plots: violin plots for nFeature_RNA, nCount_RNA, pct.mt,
#' and pct.ribo; and scatter plots of nCount vs nFeature, nCount vs pct.mt, and
#' nCount vs pct.ribo. Use the plots to determine appropriate filtering thresholds.
#'
#' @param seu A Seurat object.
#' @param mt_pattern Regex pattern to identify mitochondrial genes.
#'   Use "^mt-" for mouse, "^MT-" for human. Default "^mt-".
#' @param ribo_pattern Regex pattern to identify ribosomal genes.
#'   Use "^Rp[sl]" for mouse, "^RP[SL]" for human. Default "^Rp[sl]".
#' @param group_by Metadata column to use as the x-axis on violin plots.
#'   Default "project.name".
#'
#' @return The input Seurat object with pct.mt and pct.ribo added to metadata.
#'
#' @examples
#' seu <- plot_qc(seu, mt_pattern = "^mt-", ribo_pattern = "^Rp[sl]")
#' seu <- subset(seu, nFeature_RNA > 500 & nFeature_RNA < 6000 & pct.mt < 15)
plot_qc <- function(seu,
                    mt_pattern   = "^mt-",
                    ribo_pattern = "^Rp[sl]",
                    group_by     = "project.name") { ... }


#' @title Preprocess query and reference Seurat objects and transfer cell type labels
#'
#' @description
#' Runs a full preprocessing pipeline on both query and reference objects
#' (normalization, variable feature selection, scaling, PCA), finds transfer
#' anchors between them, and transfers the specified label from reference to query.
#' Supports both LogNormalize and SCTransform workflows.
#'
#' @param query Query Seurat object with raw counts. QC filtering should be
#'   completed before passing in.
#' @param reference Reference Seurat object with raw counts and known cell type
#'   annotations in metadata.
#' @param normalization Normalization method. One of "lognorm" (NormalizeData ->
#'   FindVariableFeatures -> ScaleData -> PCA) or "SCT" (SCTransform -> PCA).
#' @param ref_label_col Name of the metadata column in reference to transfer
#'   (e.g. "cell_subclass").
#' @param query_label_col Name of the new metadata column written into the query
#'   object containing transferred labels. Default "predicted.label".
#' @param dims Integer vector of PCA dimensions to use for FindTransferAnchors
#'   and TransferData. Default 1:30.
#' @param n_features Number of variable features to select. Default 3000.
#' @param verbose Whether to print progress from Seurat calls. Default FALSE.
#'
#' @return The query Seurat object with predicted labels in query_label_col and
#'   per-class prediction scores in prediction.score.* metadata columns.
#'
#' @examples
#' query <- preprocess_and_run_transferanchor(
#'   query           = query,
#'   reference       = ref,
#'   normalization   = "lognorm",
#'   ref_label_col   = "cell_subclass",
#'   query_label_col = "predicted_subclass",
#'   dims            = 1:30
#' )
preprocess_and_run_transferanchor <- function(query,
                                              reference,
                                              normalization   = c("lognorm", "SCT"),
                                              ref_label_col,
                                              query_label_col = "predicted.label",
                                              dims            = 1:30,
                                              n_features      = 3000,
                                              verbose         = FALSE) { ... }


#' @title Export annotated Seurat object as counts and metadata CSVs
#'
#' @description
#' Writes the raw counts matrix and cell metadata from a Seurat object to
#' separate CSV files. These files serve as the input to the downstream
#' Python analysis pipeline.
#'
#' @param seu A Seurat object.
#' @param counts_output_file Path for the output counts CSV
#'   (e.g. "data/processed/output_counts.csv").
#' @param meta_output_file Path for the output metadata CSV
#'   (e.g. "data/processed/output_metadata.csv").
#'
#' @return NULL invisibly. Files are written to the specified paths.
#'
#' @examples
#' save_annotated_data(
#'   seu                = query,
#'   counts_output_file = "data/processed/output_counts.csv",
#'   meta_output_file   = "data/processed/output_metadata.csv"
#' )
save_annotated_data <- function(seu, counts_output_file, meta_output_file) { ... }


# library(Seurat)
# library(data.table)
# library(ggplot2)
# library(patchwork)
# 
# 
# load_seurat <- function(counts_path,
#                         meta_path = NULL,
#                         transpose = FALSE,
#                         project   = "SeuratProject") {
#   
#   .read <- function(path) {
#     df <- fread(path, data.table = FALSE)
#     rownames(df) <- df[[1]]
#     df[[1]] <- NULL
#     df
#   }
#   
#   counts <- .read(counts_path)
#   if (transpose) counts <- as.data.frame(t(counts))
#   
#   meta <- NULL
#   if (!is.null(meta_path)) {
#     meta <- .read(meta_path)
#     shared <- intersect(colnames(counts), rownames(meta))
#     if (length(shared) == 0)
#       stop("No overlapping cell IDs between counts columns and metadata rownames.")
#     counts <- counts[, shared, drop = FALSE]
#     meta   <- meta[shared, , drop = FALSE]
#   }
#   
#   seu <- CreateSeuratObject(
#     counts    = counts,
#     meta.data = meta,
#     project   = project
#   )
#   message(sprintf("Created Seurat object: %d genes × %d cells", nrow(seu), ncol(seu)))
#   return(seu)
# }
# 
# plot_qc <- function(seu,
#                     mt_pattern     = "^mt-",
#                     ribo_pattern = "^Rp[sl]",
#                     group_by     = "project.name"
#                     ) {
#   seu$project.name <- seu@project.name
# 
#   seu[["pct.mt"]] <- PercentageFeatureSet(seu, pattern = mt_pattern)
#   seu[["pct.ribo"]] <- PercentageFeatureSet(seu, pattern = ribo_pattern)
#   
#   meta = seu@meta.data
#   
#   v1 <- ggplot(meta, aes(x = .data[[group_by]], y = nFeature_RNA)) + geom_violin()
#   v2 <- ggplot(meta, aes(x = .data[[group_by]], y = nCount_RNA)) + geom_violin()
#   v3 <- ggplot(meta, aes(x = .data[[group_by]], y = pct.mt)) + geom_violin()
#   v4 <- ggplot(meta, aes(x = .data[[group_by]], y = pct.ribo)) + geom_violin()
#   
#   s1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
#     NoLegend()
#   s2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "pct.mt") +
#     NoLegend()
#   s3 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "pct.ribo") + NoLegend()
#   
#   print(v1 | v2 | v3 | v4)  
#   print(s1 | s2 | s3) 
#   
#   return(seu)
# }
# 
# 
# preprocess_and_run_transferanchor <- function(query,
#                                               reference,
#                                               normalization   = c("lognorm", "SCT"),
#                                               ref_label_col,
#                                               query_label_col = "predicted.label",
#                                               dims            = 1:30,
#                                               n_features      = 3000,
#                                               verbose         = FALSE) {
#   
#   normalization <- match.arg(normalization)
#   
#   .preprocess <- function(seu) {
#     if (normalization == "lognorm") {
#       seu <- NormalizeData(seu, verbose = verbose)
#       seu <- FindVariableFeatures(seu, nfeatures = n_features, verbose = verbose)
#       seu <- ScaleData(seu, verbose = verbose)
#       seu <- RunPCA(seu, verbose = verbose)
#     } else {
#       # SCTransform: handles norm + variable features + scaling in one call
#       seu <- SCTransform(seu, variable.features.n = n_features, verbose = verbose)
#       seu <- RunPCA(seu, verbose = verbose)
#     }
#     seu
#   }
#   
#   message("Preprocessing reference...")
#   reference <- .preprocess(reference)
#   
#   message("Preprocessing query...")
#   query <- .preprocess(query)
#   
#   # ── Find anchors ────────────────────────────────────────────────────────────
#   message("Finding transfer anchors...")
#   norm_method <- if (normalization == "SCT") "SCT" else "LogNormalize"
#   
#   anchors <- FindTransferAnchors(
#     reference        = reference,
#     query            = query,
#     normalization.method = norm_method,
#     dims             = dims,
#     verbose          = verbose
#   )
#   
#   # ── Transfer labels ─────────────────────────────────────────────────────────
#   message(sprintf("Transferring '%s' → '%s'...", ref_label_col, query_label_col))
#   
#   predictions <- TransferData(
#     anchorset  = anchors,
#     refdata    = reference[[ref_label_col, drop = TRUE]],
#     dims       = dims,
#     verbose    = verbose
#   )
#   
#   # predicted.id  → user-specified column name
#   # prediction.score.* columns kept as-is
#   query <- AddMetaData(query, metadata = predictions)
#   query[[query_label_col]] <- query[["predicted.id"]]
#   
#   message("Done. Predicted label distribution:")
#   print(table(query[[query_label_col, drop = TRUE]]))
#   
#   return(query)
# }
# 
# 
# 
# save_annotated_data <- function(seu, counts_output_file, meta_output_file){
#   write.csv(as.matrix(GetAssayData(seu, slot = "counts")), 
#             counts_output_file)
#   # metadata
#   write.csv(seu@meta.data, meta_output_file)
# }


