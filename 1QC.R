###############################################################################
## Step 1: snRNA-seq filtering workflow
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

mt_pattern <- "^soybean-mt-"
cp_pattern <- "^soybean-cp-"

###############################################################################
## 3. Read 10X matrices and create Seurat objects
###############################################################################

scn.all <- vector("list", length(samples))
names(scn.all) <- samples

for (i in seq_along(samples)) {
  
  message("Reading sample: ", samples[i])
  
  counts <- Read10X(
    data.dir = sample_dirs[i]
  )
  
  scn.all[[i]] <- CreateSeuratObject(
    counts = counts,
    min.cells = 5,
    min.features = 300,
    project = samples[i]
  )
  
  scn.all[[i]]$library <- samples[i]
}

raw_cell_number <- sapply(
  scn.all,
  ncol
)

###############################################################################
## 4. Inspect raw nCount_RNA and nFeature_RNA before filtering
###############################################################################

scn_combine_raw <- merge(
  x = scn.all[[1]],
  y = scn.all[-1],
  add.cell.ids = samples
)

p <- VlnPlot(
  scn_combine_raw,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "orig.ident",
  ncol = 2,
  pt.size = 0
)

ggsave(
  filename = file.path(output_dir, "scn_combine_before_cell_filter.pdf"),
  plot = p,
  width = 20,
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
      nFeature_RNA < 10000
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

p <- VlnPlot(
  scn_combine_basic_filter,
  features = c("nFeature_RNA", "nCount_RNA"),
  group.by = "orig.ident",
  ncol = 2,
  pt.size = 0
)

ggsave(
  filename = file.path(output_dir, "scn_combine_after_nCount_nFeature_filter.pdf"),
  plot = p,
  width = 20,
  height = 6
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
  
  obj <- ScaleData(obj)
  
  obj <- RunPCA(obj)
  
  obj <- FindNeighbors(
    obj,
    dims = 1:20
  )
  
  obj <- FindClusters(
    obj,
    resolution = 1,
    verbose = FALSE
  )
  
  obj <- RunUMAP(
    obj,
    dims = 1:20
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
      bcmvn$pK[which.max(bcmvn$MeanBC)]
    )
  )
  
  annotations <- obj$seurat_clusters
  
  homotypic.prop <- modelHomotypic(
    annotations
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
  
  doublet_table <- as.data.frame(
    table(obj$doublet_info)
  )
  
  colnames(doublet_table) <- c("doublet_info", "n_cells")
  
  n_singlet <- doublet_table$n_cells[
    doublet_table$doublet_info == "Singlet"
  ]
  
  n_doublet <- doublet_table$n_cells[
    doublet_table$doublet_info == "Doublet"
  ]
  
  if (length(n_singlet) == 0) {
    n_singlet <- 0
  }
  
  if (length(n_doublet) == 0) {
    n_doublet <- 0
  }
  
  doublet_filter_summary[[i]] <- data.frame(
    sample = samples[i],
    cells_before_doublet_filter = ncol(obj),
    pK = pK_use,
    expected_doublet_rate = 0.075,
    nExp_poi = nExp_poi,
    homotypic_prop = homotypic.prop,
    nExp_poi_adjusted = nExp_poi.adj,
    n_singlet = n_singlet,
    n_doublet = n_doublet
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
  
  if (length(remove_cols) > 0) {
    
    obj@meta.data <- obj@meta.data[
      ,
      !colnames(obj@meta.data) %in% remove_cols,
      drop = FALSE
    ]
  }
  
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
## 7. Calculate mitochondrial and chloroplast proportions
###############################################################################

organelle_qc_after_doublet <- vector("list", length(samples))

for (i in seq_along(samples)) {
  
  message("Calculating organelle proportions for sample: ", samples[i])
  
  mt_features <- rownames(scn.all.filterdouble[[i]])[
    grepl(
      mt_pattern,
      rownames(scn.all.filterdouble[[i]]),
      ignore.case = TRUE
    )
  ]
  
  cp_features <- rownames(scn.all.filterdouble[[i]])[
    grepl(
      cp_pattern,
      rownames(scn.all.filterdouble[[i]]),
      ignore.case = TRUE
    )
  ]
  
  mt_features <- intersect(
    mt_features,
    rownames(scn.all.filterdouble[[i]])
  )
  
  cp_features <- intersect(
    cp_features,
    rownames(scn.all.filterdouble[[i]])
  )
  
  if (length(mt_features) > 0) {
    
    scn.all.filterdouble[[i]][["percent.mt"]] <- PercentageFeatureSet(
      scn.all.filterdouble[[i]],
      features = mt_features
    )
    
  } else {
    
    scn.all.filterdouble[[i]]$percent.mt <- 0
  }
  
  if (length(cp_features) > 0) {
    
    scn.all.filterdouble[[i]][["percent.cp"]] <- PercentageFeatureSet(
      scn.all.filterdouble[[i]],
      features = cp_features
    )
    
  } else {
    
    scn.all.filterdouble[[i]]$percent.cp <- 0
  }
  
  scn.all.filterdouble[[i]]$percent.organelle <-
    scn.all.filterdouble[[i]]$percent.mt +
    scn.all.filterdouble[[i]]$percent.cp
  
  organelle_qc_after_doublet[[i]] <- data.frame(
    sample = samples[i],
    n_cells_after_doublet = ncol(scn.all.filterdouble[[i]]),
    n_mt_genes = length(mt_features),
    n_cp_genes = length(cp_features),
    median_percent_mt = median(scn.all.filterdouble[[i]]$percent.mt),
    median_percent_cp = median(scn.all.filterdouble[[i]]$percent.cp),
    median_percent_organelle = median(scn.all.filterdouble[[i]]$percent.organelle),
    max_percent_organelle = max(scn.all.filterdouble[[i]]$percent.organelle)
  )
}

organelle_qc_after_doublet <- bind_rows(
  organelle_qc_after_doublet
)

write.csv(
  organelle_qc_after_doublet,
  file = file.path(output_dir, "organelle_percent_summary_after_doublet_filter.csv"),
  row.names = FALSE
)

scn_combine_after_doublet <- merge(
  x = scn.all.filterdouble[[1]],
  y = scn.all.filterdouble[-1],
  add.cell.ids = samples
)

p <- VlnPlot(
  scn_combine_after_doublet,
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
  filename = file.path(output_dir, "qc_after_nCount_nFeature_and_doublet_filter.pdf"),
  plot = p,
  width = 30,
  height = 6
)

rm(scn_combine_after_doublet)

###############################################################################
## 8. Classify soybean and SCN nuclei, then remove organelle genes
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
  
  colnames(species_table) <- c("species", "n_cells")
  species_table$sample <- samples[i]
  
  species_classification_summary[[i]] <- species_table[
    ,
    c("sample", "species", "n_cells")
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
    grepl("^ann1\\.Glyma", all_genes) &
      !all_genes %in% organelle_genes
  ]
  
  scn_genes <- all_genes[
    grepl("^Hetgly", all_genes) &
      !all_genes %in% organelle_genes
  ]
  
  genes_keep <- unique(
    c(
      soybean_genes,
      scn_genes
    )
  )
  
  cells_keep <- colnames(obj)[
    obj$species %in% c("soybean", "scn")
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
## 9. Merge all filtered samples
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
## 10. Export compact filtering summary
###############################################################################

filtering_cell_number_summary <- data.frame(
  sample = samples,
  raw_cells = as.integer(raw_cell_number),
  cells_after_nCount_nFeature_filter = as.integer(cell_number_after_basic_filter),
  cells_after_doublet_filter = as.integer(cell_number_after_doublet),
  cells_after_species_filter = as.integer(cell_number_after_species_filter)
)

write.csv(
  filtering_cell_number_summary,
  file = file.path(output_dir, "filtering_cell_number_summary.csv"),
  row.names = FALSE
)

###############################################################################
## 11. Inspect final QC metrics
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
## 12. Normalize, identify variable features, scale data, and run PCA
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
## 13. Save filtered objects
###############################################################################

saveRDS(
  scn_combine_filter_species,
  file = file.path(output_dir, "scn_combine_filter_species_list.rds")
)

saveRDS(
  scn_combine_filter,
  file = file.path(output_dir, "scn_combine_filter.rds")
)