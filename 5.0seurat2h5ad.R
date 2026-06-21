###############################################################################
## Step 5: Prepare transition objects and export soybean count matrix
##         for Procambium-to-Syncytium analysis
###############################################################################

rm(list = ls())
set.seed(1)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(dplyr)
library(Matrix)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "3_annotation/scn_combine_filter_harmony_celltype.rds"
output_dir <- "5_transition_objects"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load annotated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"

scn_combine_filter_harmony <- JoinLayers(
  scn_combine_filter_harmony
)

required_meta <- c(
  "seurat_clusters",
  "cell_type",
  "orig.ident",
  "breed"
)

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

scn_combine_filter_harmony$breed <- recode(
  scn_combine_filter_harmony$breed,
  "Wiliams82" = "Williams82"
)

###############################################################################
## 4. Add infection time metadata
###############################################################################

scn_combine_filter_harmony$time <- case_when(
  grepl("^[FPW]1_", scn_combine_filter_harmony$orig.ident) ~ "1 dpi",
  grepl("^[FPW]3_", scn_combine_filter_harmony$orig.ident) ~ "3 dpi",
  grepl("^[FPW]7_", scn_combine_filter_harmony$orig.ident) ~ "7 dpi",
  TRUE ~ NA_character_
)

###############################################################################
## 5. Export soybean count matrix and metadata as CSV files
###############################################################################

counts_mat <- GetAssayData(
  scn_combine_filter_harmony,
  assay = "RNA",
  layer = "counts"
)

soybean_genes <- rownames(counts_mat)[
  grepl("^(ann1\\.)?Glyma\\.", rownames(counts_mat))
]

counts_mat_soybean <- counts_mat[
  soybean_genes,
  ,
  drop = FALSE
]

write.csv(
  as.matrix(counts_mat_soybean),
  file = file.path(output_dir, "scn_combine_filter_harmony_soybean_count.csv"),
  quote = FALSE
)

metadata_soybean <- scn_combine_filter_harmony@meta.data[
  colnames(counts_mat_soybean),
  ,
  drop = FALSE
]

metadata_soybean$barcode <- rownames(metadata_soybean)

write.csv(
  metadata_soybean,
  file = file.path(output_dir, "scn_combine_filter_harmony_soybean_metadata.csv"),
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 6. Extract 1 dpi vascular and syncytium-related cell types
###############################################################################

syn_related_celltypes <- c(
  "Pericycle",
  "Procambium",
  "Phloem",
  "Xylem",
  "Syncytium"
)

syn_related_1dpi <- subset(
  scn_combine_filter_harmony,
  subset = cell_type %in% syn_related_celltypes & time == "1 dpi"
)

write.csv(
  as.data.frame(
    table(
      syn_related_1dpi$breed,
      syn_related_1dpi$cell_type
    )
  ),
  file = file.path(output_dir, "cell_number_by_breed_and_celltype_1dpi.csv"),
  row.names = FALSE
)

###############################################################################
## 7. Split and export transition objects by cultivar
###############################################################################

cultivar_objects <- SplitObject(
  syn_related_1dpi,
  split.by = "breed"
)

cultivar_objects <- cultivar_objects[
  c("Forrest", "PI88788", "Williams82")
]

for (cultivar in names(cultivar_objects)) {
  
  message("Exporting transition object for: ", cultivar)
  
  cultivar_dir <- file.path(
    output_dir,
    cultivar
  )
  
  dir.create(
    cultivar_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  obj <- cultivar_objects[[cultivar]]
  
  obj@meta.data <- obj@meta.data[
    ,
    !colnames(obj@meta.data) %in% c("doublet_info", "species"),
    drop = FALSE
  ]
  
  obj <- JoinLayers(
    obj
  )
  
  obj$barcode <- colnames(obj)
  
  if ("umap" %in% names(obj@reductions)) {
    obj$UMAP_1 <- obj@reductions$umap@cell.embeddings[, 1]
    obj$UMAP_2 <- obj@reductions$umap@cell.embeddings[, 2]
  }
  
  mat <- GetAssayData(
    obj,
    assay = "RNA",
    layer = "counts"
  )
  
  mat <- mat[
    grepl("^(ann1\\.)?Glyma\\.", rownames(mat)),
    ,
    drop = FALSE
  ]
  
  saveRDS(
    obj,
    file = file.path(
      cultivar_dir,
      paste0(cultivar, "_vascular_syn_1dpi.rds")
    )
  )
  
  Matrix::writeMM(
    mat,
    file = file.path(
      cultivar_dir,
      paste0(cultivar, "_vascular_syn_1dpi.mtx")
    )
  )
  
  write.csv(
    obj@meta.data,
    file = file.path(
      cultivar_dir,
      paste0(cultivar, "_vascular_syn_1dpi_metadata.csv")
    ),
    quote = FALSE,
    row.names = FALSE
  )
  
  gene_names <- data.frame(
    gene = sub("^ann1\\.", "", rownames(mat))
  )
  
  write.table(
    gene_names,
    file = file.path(
      cultivar_dir,
      paste0(cultivar, "_vascular_syn_1dpi_gene_names.csv")
    ),
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  if ("pca" %in% names(obj@reductions)) {
    write.csv(
      obj@reductions$pca@cell.embeddings,
      file = file.path(
        cultivar_dir,
        paste0(cultivar, "_vascular_syn_1dpi_pca.csv")
      ),
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  if ("umap" %in% names(obj@reductions)) {
    write.csv(
      obj@reductions$umap@cell.embeddings,
      file = file.path(
        cultivar_dir,
        paste0(cultivar, "_vascular_syn_1dpi_umap.csv")
      ),
      quote = FALSE,
      row.names = FALSE
    )
  }
}