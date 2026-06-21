###############################################################################
## Step 4: Pseudo-bulk differential expression analysis for Syncytium
###############################################################################

rm(list = ls())
set.seed(123)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(DESeq2)
library(dplyr)
library(Matrix)
library(Matrix.utils)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "3_annotation/scn_combine_filter_harmony_celltype.rds"
output_dir <- "4_bulk_DEG_syncytium"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load annotated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"
scn_combine_filter_harmony <- JoinLayers(scn_combine_filter_harmony)

###############################################################################
## 4. Check required metadata
###############################################################################

required_meta <- c("cell_type", "library")

missing_meta <- setdiff(
  required_meta,
  colnames(scn_combine_filter_harmony@meta.data)
)

if (length(missing_meta) > 0) {
  stop(
    "The following metadata columns are missing: ",
    paste(missing_meta, collapse = ", ")
  )
}

###############################################################################
## 5. Subset Syncytium nuclei
###############################################################################

syncytium_obj <- subset(
  scn_combine_filter_harmony,
  subset = cell_type == "Syncytium"
)

syncytium_obj <- JoinLayers(syncytium_obj)

###############################################################################
## 6. Generate pseudo-bulk count matrix by library
###############################################################################

counts <- LayerData(
  syncytium_obj,
  assay = "RNA",
  layer = "counts"
)

metadata <- syncytium_obj@meta.data[
  colnames(counts),
  ,
  drop = FALSE
]

pseudo_bulk <- Matrix.utils::aggregate.Matrix(
  t(counts),
  groupings = metadata$library,
  fun = "sum"
)

pseudo_bulk <- t(as.matrix(pseudo_bulk))
storage.mode(pseudo_bulk) <- "integer"

write.csv(
  as.data.frame(colSums(pseudo_bulk)),
  file = file.path(output_dir, "Syncytium_pseudobulk_library_size.csv")
)

###############################################################################
## 7. Define sample information and comparisons
###############################################################################

sample_info <- data.frame(
  sample_id = c(
    "F1_1", "F1_2",
    "F3_1", "F3_2",
    "F7_1", "F7_2",
    "P1_1", "P1_2",
    "P3_1", "P3_2",
    "P7_1", "P7_2",
    "W1_1", "W1_2",
    "W3_1", "W3_2",
    "W7_1", "W7_2"
  ),
  cultivar = c(
    rep("Forrest", 6),
    rep("PI88788", 6),
    rep("Williams82", 6)
  ),
  timepoint = rep(
    c("1dpi", "1dpi", "3dpi", "3dpi", "7dpi", "7dpi"),
    3
  ),
  replicate = rep(
    c("rep1", "rep2"),
    9
  ),
  stringsAsFactors = FALSE
)

comparisons <- data.frame(
  comparison_name = c(
    "Forrest_vs_Williams82_1dpi",
    "Forrest_vs_Williams82_3dpi",
    "Forrest_vs_Williams82_7dpi",
    "PI88788_vs_Williams82_1dpi",
    "PI88788_vs_Williams82_3dpi",
    "PI88788_vs_Williams82_7dpi"
  ),
  test_cultivar = c(
    "Forrest",
    "Forrest",
    "Forrest",
    "PI88788",
    "PI88788",
    "PI88788"
  ),
  ref_cultivar = rep(
    "Williams82",
    6
  ),
  timepoint = c(
    "1dpi",
    "3dpi",
    "7dpi",
    "1dpi",
    "3dpi",
    "7dpi"
  ),
  stringsAsFactors = FALSE
)

###############################################################################
## 8. Define DESeq2 comparison function
###############################################################################

run_deseq2_comparison <- function(
    count_matrix,
    sample_info,
    comparison_name,
    test_cultivar,
    ref_cultivar,
    timepoint,
    output_dir,
    min_count = 1,
    padj_cutoff = 0.05,
    log2fc_cutoff = 1
) {
  
  selected_samples <- sample_info %>%
    filter(
      timepoint == !!timepoint,
      cultivar %in% c(test_cultivar, ref_cultivar)
    ) %>%
    pull(sample_id)
  
  if (!all(selected_samples %in% colnames(count_matrix))) {
    stop(
      "Some selected samples were not found in pseudo-bulk matrix for ",
      comparison_name
    )
  }
  
  count_sub <- count_matrix[
    ,
    selected_samples,
    drop = FALSE
  ]
  
  count_sub <- count_sub[
    rowSums(count_sub) > min_count,
    ,
    drop = FALSE
  ]
  
  coldata <- sample_info %>%
    filter(sample_id %in% selected_samples) %>%
    arrange(match(sample_id, selected_samples))
  
  rownames(coldata) <- coldata$sample_id
  
  coldata$group <- factor(
    coldata$cultivar,
    levels = c(ref_cultivar, test_cultivar)
  )
  
  count_sub <- count_sub[
    ,
    rownames(coldata),
    drop = FALSE
  ]
  
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_sub),
    colData = coldata,
    design = ~ group
  )
  
  dds <- DESeq(dds)
  
  res <- results(
    dds,
    contrast = c("group", test_cultivar, ref_cultivar)
  )
  
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  
  normalized_counts <- as.data.frame(
    counts(dds, normalized = TRUE)
  )
  
  normalized_counts$gene_id <- rownames(normalized_counts)
  
  res_df <- res_df %>%
    left_join(
      normalized_counts,
      by = "gene_id"
    ) %>%
    relocate(gene_id)
  
  res_df$change <- ifelse(
    !is.na(res_df$padj) &
      res_df$padj < padj_cutoff &
      abs(res_df$log2FoldChange) > log2fc_cutoff,
    ifelse(
      res_df$log2FoldChange > log2fc_cutoff,
      "Up",
      "Down"
    ),
    "NoDiff"
  )
  
  sig_res <- res_df %>%
    filter(change %in% c("Up", "Down")) %>%
    arrange(padj)
  
  write.table(
    res_df,
    file = file.path(
      output_dir,
      paste0("full_DESeq2_Syncytium_", comparison_name, ".txt")
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  write.table(
    sig_res,
    file = file.path(
      output_dir,
      paste0("sig_DESeq2_Syncytium_", comparison_name, ".txt")
    ),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  summary_df <- data.frame(
    comparison = comparison_name,
    test_cultivar = test_cultivar,
    ref_cultivar = ref_cultivar,
    timepoint = timepoint,
    n_genes_tested = nrow(res_df),
    n_sig_genes = nrow(sig_res),
    n_up_in_test = sum(sig_res$change == "Up"),
    n_down_in_test = sum(sig_res$change == "Down")
  )
  
  return(summary_df)
}

###############################################################################
## 9. Run DESeq2 for Syncytium comparisons
###############################################################################

deg_summary <- list()

for (i in seq_len(nrow(comparisons))) {
  
  message("Running comparison: ", comparisons$comparison_name[i])
  
  deg_summary[[i]] <- run_deseq2_comparison(
    count_matrix = pseudo_bulk,
    sample_info = sample_info,
    comparison_name = comparisons$comparison_name[i],
    test_cultivar = comparisons$test_cultivar[i],
    ref_cultivar = comparisons$ref_cultivar[i],
    timepoint = comparisons$timepoint[i],
    output_dir = output_dir,
    min_count = 1,
    padj_cutoff = 0.05,
    log2fc_cutoff = 1
  )
}

deg_summary <- bind_rows(
  deg_summary
)

write.csv(
  deg_summary,
  file = file.path(output_dir, "Syncytium_DESeq2_summary.csv"),
  row.names = FALSE
)

print(deg_summary)