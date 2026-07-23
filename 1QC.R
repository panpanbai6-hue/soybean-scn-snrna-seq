###############################################################################
##snRNA-seq filtering workflow
###############################################################################

rm(list = ls())
set.seed(1234)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(DoubletFinder)

###############################################################################
## 2. Set paths and sample names
###############################################################################

matrix_dir <- "matrix_two_v4"
output_dir <- "1_filter"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "F1_1", "F1_2",
  "F3_1", "F3_2",
  "F7_1", "F7_2",
  "P1_1", "P1_2",
  "P3_1", "P3_2",
  "P7_1", "P7_2",
  "W1_1", "W1_2",
  "W3_1", "W3_2",
  "W7_1", "W7_2"
)

sample_dirs <- file.path(matrix_dir, samples)

mt_pattern <- "soybean-mt-"
cp_pattern <- "soybean-cp-"

###############################################################################
## 3. Read 10X matrices, create Seurat objects,
##    and calculate mitochondrial and chloroplast proportions
###############################################################################

scn.all <- vector("list", length(samples))
names(scn.all) <- samples

for (i in seq_along(samples)) {
  
  message("Reading sample: ", samples[i])
  
  counts <- Read10X(
    data.dir = sample_dirs[i]
  )
  
  obj <- CreateSeuratObject(
    counts = counts,
    min.cells = 5,
    min.features = 300,
    project = samples[i]
  )
  
  obj$library <- samples[i]
  
  obj[["percent.mt"]] <- PercentageFeatureSet(
    obj,
    pattern = mt_pattern
  )
  
  obj[["percent.cp"]] <- PercentageFeatureSet(
    obj,
    pattern = cp_pattern
  )
  
  obj$percent.organelle <- obj$percent.mt + obj$percent.cp
  
  scn.all[[i]] <- obj
}

raw_cell_number <- sapply(
  scn.all,
  ncol
)

###############################################################################
## 4. Inspect raw nCount_RNA, nFeature_RNA,
##    mitochondrial and chloroplast proportions before filtering
###############################################################################

scn_combine_raw <- merge(
  x = scn.all[[1]],
  y = scn.all[-1],
  add.cell.ids = samples
)

p <- VlnPlot(
  scn_combine_raw,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt",
    "percent.cp",
    "percent.organelle"
  ),
  group.by = "orig.ident",
  ncol = 5,
  pt.size = 0
)

ggsave(
  filename = file.path(output_dir, "scn_combine_before_cell_filter_qc.pdf"),
  plot = p,
  width = 30,
  height = 6
)

rm(scn_combine_raw)

###############################################################################
## 5. Filter nuclei by nCount_RNA and nFeature_RNA
###############################################################################

scn.all.filter <- vector("list", length(samples))
names(scn.all.filter) <- samples

## F1, F3, F7, P1, P3, W1 and W3 samples.
for (i in c(1:10, 13:16)) {
  
  scn.all.filter[[i]] <- subset(
    scn.all[[i]],
    subset = nCount_RNA < 25000 &
      nFeature_RNA > 450 &
      nFeature_RNA < 11000
  )
}

## P7 and W7 samples.
for (i in c(11, 12, 17, 18)) {
  
  scn.all.filter[[i]] <- subset(
    scn.all[[i]],
    subset = nCount_RNA < 15000 &
      nFeature_RNA > 450 &
      nFeature_RNA < 7000
  )
}

cell_number_after_basic_filter <- sapply(
  scn.all.filter,
  ncol
)

scn_combine_basic_filter <- merge(
  x = scn.all.filter[[1]],
  y = scn.all.filter[-1],
  add.cell.ids = samples
)

rm(scn_combine_basic_filter)

###############################################################################
## 6. Remove doublets using DoubletFinder
###############################################################################

scn.all.filterdouble <- vector("list", length(samples))
names(scn.all.filterdouble) <- samples

doublet_filter_summary <- vector("list", length(samples))

for (i in seq_along(samples)) {
  
  message("Running DoubletFinder for sample: ", samples[i])
  
  set.seed(1000 + i)
  
  obj <- scn.all.filter[[i]]
  
  obj <- NormalizeData(
    obj,
    verbose = FALSE
  )
  
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  obj <- ScaleData(
    obj,
    verbose = FALSE
  )
  
  obj <- RunPCA(
    obj,
    verbose = FALSE
  )
  
  obj <- FindNeighbors(
    obj,
    dims = 1:20,
    verbose = FALSE
  )
  
  obj <- FindClusters(
    obj,
    resolution = 1,
    verbose = FALSE
  )
  
  obj <- RunUMAP(
    obj,
    dims = 1:20,
    verbose = FALSE
  )
  
  sweep.res_list <- paramSweep(
    obj,
    PCs = 1:20,
    sct = FALSE
  )
  
  sweep.stats <- summarizeSweep(
    sweep.res_list,
    GT = FALSE
  )
  
  bcmvn <- find.pK(
    sweep.stats
  )
  
  pK_use <- as.numeric(
    as.character(
      bcmvn$pK[
        which.max(bcmvn$MeanBC)
      ]
    )
  )
  
  homotypic.prop <- modelHomotypic(
    obj$seurat_clusters
  )
  
  nExp_poi <- round(
    0.075 * ncol(obj)
  )
  
  nExp_poi.adj <- round(
    nExp_poi * (1 - homotypic.prop)
  )
  
  set.seed(1000 + i)
  
  obj <- doubletFinder(
    obj,
    PCs = 1:20,
    pN = 0.25,
    pK = pK_use,
    nExp = nExp_poi.adj,
    reuse.pANN = NULL,
    sct = FALSE
  )
  
  doublet_cols <- grep(
    "^DF.classifications",
    colnames(obj@meta.data),
    value = TRUE
  )
  
  if (length(doublet_cols) == 0) {
    stop(
      "DoubletFinder classification column was not found for sample: ",
      samples[i]
    )
  }
  
  obj$doublet_info <- obj@meta.data[[doublet_cols[length(doublet_cols)]]]
  
  doublet_filter_summary[[i]] <- data.frame(
    sample = samples[i],
    cells_before_doublet_filter = ncol(obj),
    pK = pK_use,
    expected_doublet_rate = 0.075,
    nExp_poi = nExp_poi,
    homotypic_prop = homotypic.prop,
    nExp_poi_adjusted = nExp_poi.adj,
    n_singlet = sum(obj$doublet_info == "Singlet"),
    n_doublet = sum(obj$doublet_info == "Doublet")
  )
  
  obj <- subset(
    obj,
    subset = doublet_info == "Singlet"
  )
  
  obj$library <- samples[i]
  obj$barcode_raw <- colnames(obj)
  
  pann_cols <- grep(
    "^pANN_",
    colnames(obj@meta.data),
    value = TRUE
  )
  
  remove_cols <- c(
    pann_cols,
    doublet_cols
  )
  
  obj@meta.data <- obj@meta.data[
    ,
    !colnames(obj@meta.data) %in% remove_cols,
    drop = FALSE
  ]
  
  scn.all.filterdouble[[i]] <- obj
}

doublet_filter_summary <- bind_rows(
  doublet_filter_summary
)

write.csv(
  doublet_filter_summary,
  file = file.path(output_dir, "doublet_filter_summary.csv"),
  row.names = FALSE
)

cell_number_after_doublet <- sapply(
  scn.all.filterdouble,
  ncol
)

###############################################################################
## 7. Classify soybean and SCN nuclei, then remove organelle genes
###############################################################################

scn_combine_filter_species <- vector("list", length(samples))
names(scn_combine_filter_species) <- samples

species_classification_summary <- vector("list", length(samples))

for (i in seq_along(samples)) {
  
  message("Classifying species and removing organelle genes for sample: ", samples[i])
  
  obj <- scn.all.filterdouble[[i]]
  
  obj[["prop_soybean"]] <- PercentageFeatureSet(
    obj,
    pattern = "ann1\\.Glyma"
  )
  
  obj[["prop_scn"]] <- PercentageFeatureSet(
    obj,
    pattern = "Hetgly"
  )
  
  obj$species <- case_when(
    obj$prop_soybean >= 80 ~ "soybean",
    obj$prop_scn >= 80 ~ "scn",
    TRUE ~ "mixed"
  )
  
  species_table <- as.data.frame(
    table(obj$species)
  )
  
  colnames(species_table) <- c(
    "species",
    "n_cells"
  )
  
  species_table$sample <- samples[i]
  
  species_classification_summary[[i]] <- species_table[
    ,
    c(
      "sample",
      "species",
      "n_cells"
    )
  ]
  
  all_genes <- rownames(obj)
  
  organelle_genes <- all_genes[
    grepl(
      mt_pattern,
      all_genes,
      ignore.case = TRUE
    ) |
      grepl(
        cp_pattern,
        all_genes,
        ignore.case = TRUE
      )
  ]
  
  soybean_genes <- all_genes[
    grepl(
      "^ann1\\.Glyma",
      all_genes
    ) &
      !all_genes %in% organelle_genes
  ]
  
  scn_genes <- all_genes[
    grepl(
      "^Hetgly",
      all_genes
    ) &
      !all_genes %in% organelle_genes
  ]
  
  genes_keep <- unique(
    c(
      soybean_genes,
      scn_genes
    )
  )
  
  cells_keep <- colnames(obj)[
    obj$species %in% c(
      "soybean",
      "scn"
    )
  ]
  
  filtered_obj <- subset(
    obj,
    cells = cells_keep,
    features = genes_keep
  )
  
  filtered_obj <- JoinLayers(
    filtered_obj
  )
  
  count_mat <- LayerData(
    filtered_obj,
    assay = "RNA",
    layer = "counts"
  )
  
  filtered_obj$nCount_RNA <- Matrix::colSums(
    count_mat
  )
  
  filtered_obj$nFeature_RNA <- Matrix::colSums(
    count_mat > 0
  )
  
  filtered_obj$library <- samples[i]
  filtered_obj$barcode_raw <- colnames(filtered_obj)
  
  scn_combine_filter_species[[i]] <- filtered_obj
}

species_classification_summary <- bind_rows(
  species_classification_summary
)

write.csv(
  species_classification_summary,
  file = file.path(output_dir, "species_classification_summary.csv"),
  row.names = FALSE
)

cell_number_after_species_filter <- sapply(
  scn_combine_filter_species,
  ncol
)

###############################################################################
## 8. Merge all filtered samples
###############################################################################

scn_combine_filter <- merge(
  x = scn_combine_filter_species[[1]],
  y = scn_combine_filter_species[-1],
  add.cell.ids = samples
)

scn_combine_filter <- JoinLayers(
  scn_combine_filter
)

count_mat <- LayerData(
  scn_combine_filter,
  assay = "RNA",
  layer = "counts"
)

scn_combine_filter$nCount_RNA <- Matrix::colSums(
  count_mat
)

scn_combine_filter$nFeature_RNA <- Matrix::colSums(
  count_mat > 0
)


###############################################################################
## 9. Inspect final QC metrics
###############################################################################

p <- VlnPlot(
  scn_combine_filter,
  features = c(
    "nFeature_RNA",
    "nCount_RNA"
  ),
  group.by = "orig.ident",
  ncol = 2,
  pt.size = 0
)

ggsave(
  filename = file.path(output_dir, "scn_combine_after_all_filtering_qc.pdf"),
  plot = p,
  width = 20,
  height = 6
)

p_soybean <- VlnPlot(
  scn_combine_filter,
  features = "prop_soybean",
  group.by = "orig.ident",
  ncol = 1,
  pt.size = 0
) +
  geom_hline(
    yintercept = 80,
    linetype = "dotted",
    color = "red"
  )

ggsave(
  filename = file.path(output_dir, "prop_soybean_after_filter.pdf"),
  plot = p_soybean,
  width = 10,
  height = 6
)

p_scn <- VlnPlot(
  scn_combine_filter,
  features = "prop_scn",
  group.by = "orig.ident",
  ncol = 1,
  pt.size = 0
) +
  geom_hline(
    yintercept = 80,
    linetype = "dotted",
    color = "#284184"
  )

ggsave(
  filename = file.path(output_dir, "prop_scn_after_filter.pdf"),
  plot = p_scn,
  width = 10,
  height = 6
)

###############################################################################
## 10. Normalize, identify variable features, scale data, and run PCA
###############################################################################

scn_combine_filter <- NormalizeData(
  scn_combine_filter
)

scn_combine_filter <- FindVariableFeatures(
  scn_combine_filter,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = FALSE
)

scn_combine_filter <- ScaleData(
  scn_combine_filter
)

scn_combine_filter <- RunPCA(
  scn_combine_filter,
  npcs = 50
)

pdf(
  file = file.path(output_dir, "ElbowPlot_plot.pdf"),
  width = 8,
  height = 6
)

print(
  ElbowPlot(
    scn_combine_filter,
    ndims = 50
  )
)

dev.off()

###############################################################################
## 11. Save filtered objects
###############################################################################

saveRDS(
  scn_combine_filter_species,
  file = file.path(output_dir, "scn_combine_filter_species_list.rds")
)

saveRDS(
  scn_combine_filter,
  file = file.path(output_dir, "scn_combine_filter.rds")
)
