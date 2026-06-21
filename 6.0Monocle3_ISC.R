###############################################################################
## Step 6: Monocle3 trajectory analysis from Procambium to Syncytium
##         and syncytium-related cluster correlation analysis
###############################################################################

rm(list = ls())
set.seed(1)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(monocle3)
library(dplyr)
library(ggplot2)
library(Matrix)
library(ClusterGVis)
library(pheatmap)
library(clusterProfiler)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "3_annotation/scn_combine_filter_harmony_celltype.rds"
output_dir <- "6_monocle3_procambium_to_syncytium"

syncytium_specific_gene_file <- "5.1cellex/syncytium_specific_genes.csv"

go_annotation_file <- "annotation/GOannotation_v4.tsv"
go_info_file <- "annotation/go.tb"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Forrest"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "PI88788"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Williams82"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "GO_enrichment"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "syncytium_related_cluster_correlation"), recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load annotated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"

scn_combine_filter_harmony <- JoinLayers(
  scn_combine_filter_harmony
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

scn_combine_filter_harmony$seurat_clusters <- as.character(
  scn_combine_filter_harmony$seurat_clusters
)

###############################################################################
## 5. Syncytium-related cluster correlation analysis
###############################################################################

correlation_dir <- file.path(
  output_dir,
  "syncytium_related_cluster_correlation"
)

syn_spe_gene <- read.csv(
  syncytium_specific_gene_file,
  stringsAsFactors = FALSE
)

syn_spe_gene <- syn_spe_gene$gene

syn_related_clusters <- c(
  "13",
  "11",
  "5", "9", "14", "24",
  "21",
  "17"
)

syn_related_1dpi <- subset(
  scn_combine_filter_harmony,
  subset = seurat_clusters %in% syn_related_clusters & time == "1 dpi"
)

write.csv(
  as.data.frame(
    table(
      syn_related_1dpi$breed,
      syn_related_1dpi$seurat_clusters
    )
  ),
  file = file.path(correlation_dir, "cell_number_by_breed_and_cluster_1dpi.csv"),
  row.names = FALSE
)

###############################################################################
## 5.1 Forrest correlation analysis
###############################################################################

Forrest_cor <- subset(
  syn_related_1dpi,
  subset = breed == "Forrest"
)

av <- AggregateExpression(
  Forrest_cor,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"
)

av <- av[[1]]
rownames(av) <- sub("^ann1\\.", "", rownames(av))

syn_spe_gene_Forrest <- intersect(
  syn_spe_gene,
  rownames(av)
)

sync_mat <- as.data.frame(
  av[syn_spe_gene_Forrest, , drop = FALSE]
)

cluster17_col <- intersect(
  c("g17", "17"),
  colnames(sync_mat)
)[1]

keep_expr <- sync_mat[, cluster17_col] > 30

sync_genes_filtered <- rownames(
  sync_mat[keep_expr, , drop = FALSE]
)

cg_data <- as.data.frame(
  av[sync_genes_filtered, , drop = FALSE]
)

write.table(
  sync_genes_filtered,
  file = file.path(
    correlation_dir,
    "Forrest_syncytium_specific_genes_cluster17_expr_gt30.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

pheatmap(
  cor(cg_data, method = "pearson"),
  color = colorRampPalette(c("#f0ece1", "#FFE4E1", "#e45b5b"))(50),
  border_color = "#c0c0c0",
  display_numbers = TRUE,
  filename = file.path(
    correlation_dir,
    "Forrest_syncytium_related_cluster_heatmap_syn_spe_1dpi.pdf"
  ),
  width = 5.5,
  height = 5
)

###############################################################################
## 5.2 PI88788 correlation analysis
###############################################################################

PI88788_cor <- subset(
  syn_related_1dpi,
  subset = breed == "PI88788"
)

av <- AggregateExpression(
  PI88788_cor,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"
)

av <- av[[1]]
rownames(av) <- sub("^ann1\\.", "", rownames(av))

syn_spe_gene_PI88788 <- intersect(
  syn_spe_gene,
  rownames(av)
)

sync_mat <- as.data.frame(
  av[syn_spe_gene_PI88788, , drop = FALSE]
)

cluster17_col <- intersect(
  c("g17", "17"),
  colnames(sync_mat)
)[1]

keep_expr <- sync_mat[, cluster17_col] > 30

sync_genes_filtered <- rownames(
  sync_mat[keep_expr, , drop = FALSE]
)

cg_data <- as.data.frame(
  av[sync_genes_filtered, , drop = FALSE]
)

write.table(
  sync_genes_filtered,
  file = file.path(
    correlation_dir,
    "PI88788_syncytium_specific_genes_cluster17_expr_gt30.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

pheatmap(
  cor(cg_data, method = "pearson"),
  color = colorRampPalette(c("#f0ece1", "#FFE4E1", "#e45b5b"))(50),
  border_color = "#c0c0c0",
  display_numbers = TRUE,
  filename = file.path(
    correlation_dir,
    "PI88788_syncytium_related_cluster_heatmap_syn_spe_1dpi.pdf"
  ),
  width = 5.5,
  height = 5
)

###############################################################################
## 5.3 Williams82 correlation analysis
###############################################################################

Williams82_cor <- subset(
  syn_related_1dpi,
  subset = breed == "Williams82"
)

av <- AggregateExpression(
  Williams82_cor,
  group.by = "seurat_clusters",
  assays = "RNA",
  slot = "data"
)

av <- av[[1]]
rownames(av) <- sub("^ann1\\.", "", rownames(av))

syn_spe_gene_Williams82 <- intersect(
  syn_spe_gene,
  rownames(av)
)

sync_mat <- as.data.frame(
  av[syn_spe_gene_Williams82, , drop = FALSE]
)

cluster17_col <- intersect(
  c("g17", "17"),
  colnames(sync_mat)
)[1]

keep_expr <- sync_mat[, cluster17_col] > 30

sync_genes_filtered <- rownames(
  sync_mat[keep_expr, , drop = FALSE]
)

cg_data <- as.data.frame(
  av[sync_genes_filtered, , drop = FALSE]
)

write.table(
  sync_genes_filtered,
  file = file.path(
    correlation_dir,
    "Williams82_syncytium_specific_genes_cluster17_expr_gt30.txt"
  ),
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

pheatmap(
  cor(cg_data, method = "pearson"),
  color = colorRampPalette(c("#f0ece1", "#FFE4E1", "#e45b5b"))(50),
  border_color = "#c0c0c0",
  display_numbers = TRUE,
  filename = file.path(
    correlation_dir,
    "Williams82_syncytium_related_cluster_heatmap_syn_spe_1dpi.pdf"
  ),
  width = 5.5,
  height = 5
)

###############################################################################
## 6. Extract Procambium and Syncytium nuclei at 1 dpi for Monocle3
###############################################################################

syn_procambium_1dpi <- subset(
  scn_combine_filter_harmony,
  subset = cell_type %in% c("Procambium", "Syncytium") & time == "1 dpi"
)

pre_pseudotime_matrix <- getFromNamespace(
  "pre_pseudotime_matrix",
  "ClusterGVis"
)

###############################################################################
## 7. Forrest Monocle3 trajectory analysis
###############################################################################

Forrest <- subset(
  syn_procambium_1dpi,
  subset = breed == "Forrest"
)

DefaultAssay(Forrest) <- "RNA"

Forrest <- JoinLayers(Forrest)

Forrest_dir <- file.path(output_dir, "Forrest")

Forrest_syn_procambium_1dpi <- subset(
  Forrest,
  subset = cell_type %in% c("Procambium", "Syncytium")
)

Forrest_syn_procambium_1dpi <- JoinLayers(
  Forrest_syn_procambium_1dpi
)

expr_matrix <- as.matrix(
  LayerData(
    Forrest_syn_procambium_1dpi,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grepl("^(ann1\\.)?Glyma\\.", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- Forrest_syn_procambium_1dpi@meta.data[
  colnames(expr_matrix),
  ,
  drop = FALSE
]

gene_annotation <- data.frame(
  gene_short_name = sub("^ann1\\.", "", rownames(expr_matrix))
)

rownames(gene_annotation) <- rownames(expr_matrix)

Forrest_cds_syn_procambium_1dpi <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

Forrest_cds_syn_procambium_1dpi <- preprocess_cds(
  Forrest_cds_syn_procambium_1dpi,
  num_dim = 10
)

Forrest_cds_syn_procambium_1dpi <- align_cds(
  Forrest_cds_syn_procambium_1dpi,
  alignment_group = "orig.ident"
)

Forrest_cds_syn_procambium_1dpi <- reduce_dimension(
  Forrest_cds_syn_procambium_1dpi,
  preprocess_method = "PCA"
)

Forrest_cds_syn_procambium_1dpi <- cluster_cells(
  Forrest_cds_syn_procambium_1dpi,
  reduction_method = "UMAP"
)

p <- plot_cells(
  Forrest_cds_syn_procambium_1dpi,
  reduction_method = "UMAP",
  cell_size = 1
)

ggsave(
  file.path(Forrest_dir, "Forrest_cds_syn_procambium_1dpi_recluster.pdf"),
  plot = p,
  height = 3,
  width = 4
)

Forrest_cds_syn_procambium_1dpi <- learn_graph(
  Forrest_cds_syn_procambium_1dpi,
  use_partition = FALSE
)

p <- plot_cells(
  Forrest_cds_syn_procambium_1dpi,
  color_cells_by = "cell_type",
  trajectory_graph_color = "#ff4e50",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "Procambium" = "#C03B4C",
      "Syncytium" = "#456F97"
    )
  )

ggsave(
  file.path(Forrest_dir, "Forrest_cds_syn_procambium_1dpi_cell_type.pdf"),
  plot = p,
  height = 3,
  width = 4
)

Forrest_cds_syn_procambium_1dpi <- order_cells(
  Forrest_cds_syn_procambium_1dpi
)

p <- plot_cells(
  Forrest_cds_syn_procambium_1dpi,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  file.path(Forrest_dir, "Forrest_cds_syn_procambium_1dpi_pseudotime.pdf"),
  plot = p,
  height = 3,
  width = 5
)

Forrest_graph_test <- graph_test(
  Forrest_cds_syn_procambium_1dpi,
  neighbor_graph = "principal_graph",
  cores = 4
)

Forrest_pr_deg_ids <- rownames(
  subset(
    Forrest_graph_test,
    q_value <= 0.01 & morans_I > 0.1
  )
)

write.table(
  Forrest_pr_deg_ids,
  file = file.path(Forrest_dir, "Forrest_syn_procambium_1dpi_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  Forrest_graph_test,
  file = file.path(Forrest_dir, "Forrest_monocle3_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

Forrest_mat <- pre_pseudotime_matrix(
  cds_obj = Forrest_cds_syn_procambium_1dpi,
  gene_list = Forrest_pr_deg_ids
)

Forrest_ck <- clusterData(
  obj = Forrest_mat,
  clusterMethod = "kmeans",
  clusterNum = 3
)

Forrest_mark_genes <- character(0)

if (file.exists("forrest_TF.txt")) {
  Forrest_tf <- read.table(
    "forrest_TF.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  Forrest_mark_genes <- as.vector(Forrest_tf$V1)
}

pdf(
  file.path(Forrest_dir, "Forrest_cds_syn_procambium_1dpi_pseudotime_heatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = Forrest_ck,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = Forrest_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(2, 3, 1)
)

dev.off()

Forrest_heatmap_genes <- as.data.frame(
  Forrest_ck$long.res
)

Forrest_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  Forrest_heatmap_genes$gene
)

write.table(
  Forrest_heatmap_genes,
  file = file.path(Forrest_dir, "Forrest_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 8. PI88788 Monocle3 trajectory analysis
###############################################################################

PI88788 <- subset(
  syn_procambium_1dpi,
  subset = breed == "PI88788"
)

DefaultAssay(PI88788) <- "RNA"

PI88788 <- JoinLayers(
  PI88788
)

PI88788_dir <- file.path(output_dir, "PI88788")

PI88788_syn_procambium_1dpi <- subset(
  PI88788,
  subset = cell_type %in% c("Procambium", "Syncytium")
)

PI88788_syn_procambium_1dpi <- JoinLayers(
  PI88788_syn_procambium_1dpi
)

expr_matrix <- as.matrix(
  LayerData(
    PI88788_syn_procambium_1dpi,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grepl("^(ann1\\.)?Glyma\\.", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- PI88788_syn_procambium_1dpi@meta.data[
  colnames(expr_matrix),
  ,
  drop = FALSE
]

gene_annotation <- data.frame(
  gene_short_name = sub("^ann1\\.", "", rownames(expr_matrix))
)

rownames(gene_annotation) <- rownames(expr_matrix)

PI88788_cds_syn_procambium_1dpi <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

PI88788_cds_syn_procambium_1dpi <- preprocess_cds(
  PI88788_cds_syn_procambium_1dpi,
  num_dim = 7
)

PI88788_cds_syn_procambium_1dpi <- align_cds(
  PI88788_cds_syn_procambium_1dpi,
  alignment_group = "orig.ident"
)

PI88788_cds_syn_procambium_1dpi <- reduce_dimension(
  PI88788_cds_syn_procambium_1dpi,
  preprocess_method = "PCA"
)

PI88788_cds_syn_procambium_1dpi <- cluster_cells(
  PI88788_cds_syn_procambium_1dpi,
  reduction_method = "UMAP"
)

p <- plot_cells(
  PI88788_cds_syn_procambium_1dpi,
  reduction_method = "UMAP",
  cell_size = 1
)

ggsave(
  file.path(PI88788_dir, "PI88788_cds_syn_procambium_1dpi_recluster.pdf"),
  plot = p,
  height = 3,
  width = 4
)

PI88788_cds_syn_procambium_1dpi <- learn_graph(
  PI88788_cds_syn_procambium_1dpi
)

p <- plot_cells(
  PI88788_cds_syn_procambium_1dpi,
  color_cells_by = "cell_type",
  trajectory_graph_color = "#ff4e50",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "Procambium" = "#C03B4C",
      "Syncytium" = "#456F97"
    )
  )

ggsave(
  file.path(PI88788_dir, "PI88788_cds_syn_procambium_1dpi_cell_type.pdf"),
  plot = p,
  height = 3,
  width = 4
)

PI88788_cds_syn_procambium_1dpi <- order_cells(
  PI88788_cds_syn_procambium_1dpi
)

p <- plot_cells(
  PI88788_cds_syn_procambium_1dpi,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  file.path(PI88788_dir, "PI88788_cds_syn_procambium_1dpi_pseudotime.pdf"),
  plot = p,
  height = 3,
  width = 5
)

PI88788_graph_test <- graph_test(
  PI88788_cds_syn_procambium_1dpi,
  neighbor_graph = "principal_graph",
  cores = 4
)

PI88788_pr_deg_ids <- rownames(
  subset(
    PI88788_graph_test,
    q_value <= 0.01 & morans_I > 0.1
  )
)

write.table(
  PI88788_pr_deg_ids,
  file = file.path(PI88788_dir, "PI88788_syn_procambium_1dpi_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  PI88788_graph_test,
  file = file.path(PI88788_dir, "PI88788_monocle3_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

PI88788_mat <- pre_pseudotime_matrix(
  cds_obj = PI88788_cds_syn_procambium_1dpi,
  gene_list = PI88788_pr_deg_ids
)

PI88788_ck <- clusterData(
  obj = PI88788_mat,
  clusterMethod = "kmeans",
  clusterNum = 3
)

PI88788_mark_genes <- character(0)

if (file.exists("PI88788_TF.txt")) {
  PI88788_tf <- read.table(
    "PI88788_TF.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  PI88788_mark_genes <- as.vector(PI88788_tf$V1)
}

pdf(
  file.path(PI88788_dir, "PI88788_cds_syn_procambium_1dpi_pseudotime_heatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = PI88788_ck,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = PI88788_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(3, 2, 1)
)

dev.off()

PI88788_heatmap_genes <- as.data.frame(
  PI88788_ck$long.res
)

PI88788_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  PI88788_heatmap_genes$gene
)

write.table(
  PI88788_heatmap_genes,
  file = file.path(PI88788_dir, "PI88788_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 9. Williams82 Monocle3 trajectory analysis
###############################################################################

Williams82 <- subset(
  syn_procambium_1dpi,
  subset = breed == "Williams82"
)

DefaultAssay(Williams82) <- "RNA"

Williams82 <- JoinLayers(
  Williams82
)

Williams82_dir <- file.path(output_dir, "Williams82")

Williams82_syn_procambium_1dpi <- subset(
  Williams82,
  subset = cell_type %in% c("Procambium", "Syncytium")
)

Williams82_syn_procambium_1dpi <- JoinLayers(
  Williams82_syn_procambium_1dpi
)

expr_matrix <- as.matrix(
  LayerData(
    Williams82_syn_procambium_1dpi,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grepl("^(ann1\\.)?Glyma\\.", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- Williams82_syn_procambium_1dpi@meta.data[
  colnames(expr_matrix),
  ,
  drop = FALSE
]

gene_annotation <- data.frame(
  gene_short_name = sub("^ann1\\.", "", rownames(expr_matrix))
)

rownames(gene_annotation) <- rownames(expr_matrix)

Williams82_cds_syn_procambium_1dpi <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

Williams82_cds_syn_procambium_1dpi <- preprocess_cds(
  Williams82_cds_syn_procambium_1dpi,
  num_dim = 8
)

Williams82_cds_syn_procambium_1dpi <- align_cds(
  Williams82_cds_syn_procambium_1dpi,
  alignment_group = "orig.ident"
)

Williams82_cds_syn_procambium_1dpi <- reduce_dimension(
  Williams82_cds_syn_procambium_1dpi,
  preprocess_method = "PCA"
)

Williams82_cds_syn_procambium_1dpi <- cluster_cells(
  Williams82_cds_syn_procambium_1dpi,
  reduction_method = "UMAP"
)

p <- plot_cells(
  Williams82_cds_syn_procambium_1dpi,
  reduction_method = "UMAP",
  cell_size = 1
)

ggsave(
  file.path(Williams82_dir, "Williams82_cds_syn_procambium_1dpi_recluster.pdf"),
  plot = p,
  height = 3,
  width = 4
)

Williams82_cds_syn_procambium_1dpi <- learn_graph(
  Williams82_cds_syn_procambium_1dpi
)

p <- plot_cells(
  Williams82_cds_syn_procambium_1dpi,
  color_cells_by = "cell_type",
  trajectory_graph_color = "#ff4e50",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = TRUE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "Procambium" = "#C03B4C",
      "Syncytium" = "#456F97"
    )
  )

ggsave(
  file.path(Williams82_dir, "Williams82_cds_syn_procambium_1dpi_cell_type.pdf"),
  plot = p,
  height = 3,
  width = 4
)

Williams82_cds_syn_procambium_1dpi <- order_cells(
  Williams82_cds_syn_procambium_1dpi
)

p <- plot_cells(
  Williams82_cds_syn_procambium_1dpi,
  color_cells_by = "pseudotime",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  file.path(Williams82_dir, "Williams82_cds_syn_procambium_1dpi_pseudotime.pdf"),
  plot = p,
  height = 3,
  width = 5
)

Williams82_graph_test <- graph_test(
  Williams82_cds_syn_procambium_1dpi,
  neighbor_graph = "principal_graph",
  cores = 4
)

Williams82_pr_deg_ids <- rownames(
  subset(
    Williams82_graph_test,
    q_value <= 0.01 & morans_I > 0.1
  )
)

write.table(
  Williams82_pr_deg_ids,
  file = file.path(Williams82_dir, "Williams82_syn_procambium_1dpi_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  Williams82_graph_test,
  file = file.path(Williams82_dir, "Williams82_monocle3_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

Williams82_mat <- pre_pseudotime_matrix(
  cds_obj = Williams82_cds_syn_procambium_1dpi,
  gene_list = Williams82_pr_deg_ids
)

Williams82_ck <- clusterData(
  obj = Williams82_mat,
  clusterMethod = "kmeans",
  clusterNum = 3
)

Williams82_mark_genes <- character(0)

if (file.exists("W82_tf.txt")) {
  W82_tf <- read.table(
    "W82_tf.txt",
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  Williams82_mark_genes <- as.vector(W82_tf$V1)
}

pdf(
  file.path(Williams82_dir, "Williams82_cds_syn_procambium_1dpi_pseudotime_heatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = Williams82_ck,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = Williams82_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(3, 1, 2)
)

dev.off()

Williams82_heatmap_genes <- as.data.frame(
  Williams82_ck$long.res
)

Williams82_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  Williams82_heatmap_genes$gene
)

write.table(
  Williams82_heatmap_genes,
  file = file.path(Williams82_dir, "Williams82_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


###############################################################################
## 10. GO enrichment for pseudotime gene modules
###############################################################################

GOannotation <- read.delim(
  go_annotation_file,
  stringsAsFactors = FALSE,
  header = TRUE
)

GOinfo <- read.delim(
  go_info_file,
  stringsAsFactors = FALSE,
  header = FALSE
)

###############################################################################
## 10.1 Forrest GO enrichment
###############################################################################

Forrest_heatmap_gene_file <- file.path(
  Forrest_dir,
  "Forrest_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"
)

Forrest_heatmap_gene_table <- read.table(
  Forrest_heatmap_gene_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

for (module_id in sort(unique(Forrest_heatmap_gene_table$cluster))) {
  
  marker <- unique(
    Forrest_heatmap_gene_table$gene_clean[
      Forrest_heatmap_gene_table$cluster == module_id
    ]
  )
  
  BP <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "biological_process",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  CC <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "cellular_component",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  MF <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "molecular_function",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  BP <- as.data.frame(BP) %>%
    mutate(Category = "BP")
  
  CC <- as.data.frame(CC) %>%
    mutate(Category = "CC")
  
  MF <- as.data.frame(MF) %>%
    mutate(Category = "MF")
  
  all_go <- bind_rows(BP, CC, MF)
  
  write.table(
    all_go,
    file = file.path(
      output_dir,
      "GO_enrichment",
      paste0("Forrest_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 10.2 PI88788 GO enrichment
###############################################################################

PI88788_heatmap_gene_file <- file.path(
  PI88788_dir,
  "PI88788_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"
)

PI88788_heatmap_gene_table <- read.table(
  PI88788_heatmap_gene_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

for (module_id in sort(unique(PI88788_heatmap_gene_table$cluster))) {
  
  marker <- unique(
    PI88788_heatmap_gene_table$gene_clean[
      PI88788_heatmap_gene_table$cluster == module_id
    ]
  )
  
  BP <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "biological_process",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  CC <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "cellular_component",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  MF <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "molecular_function",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  BP <- as.data.frame(BP) %>%
    mutate(Category = "BP")
  
  CC <- as.data.frame(CC) %>%
    mutate(Category = "CC")
  
  MF <- as.data.frame(MF) %>%
    mutate(Category = "MF")
  
  all_go <- bind_rows(BP, CC, MF)
  
  write.table(
    all_go,
    file = file.path(
      output_dir,
      "GO_enrichment",
      paste0("PI88788_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 10.3 Williams82 GO enrichment
###############################################################################

Williams82_heatmap_gene_file <- file.path(
  Williams82_dir,
  "Williams82_cds_syn_procambium_1dpi_pseudotime_heatmap_genes.txt"
)

Williams82_heatmap_gene_table <- read.table(
  Williams82_heatmap_gene_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

for (module_id in sort(unique(Williams82_heatmap_gene_table$cluster))) {
  
  marker <- unique(
    Williams82_heatmap_gene_table$gene_clean[
      Williams82_heatmap_gene_table$cluster == module_id
    ]
  )
  
  BP <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "biological_process",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  CC <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "cellular_component",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  MF <- enricher(
    marker,
    TERM2GENE = GOannotation[
      GOannotation$level == "molecular_function",
      c(2, 1)
    ],
    TERM2NAME = GOinfo[, 1:2],
    qvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  BP <- as.data.frame(BP) %>%
    mutate(Category = "BP")
  
  CC <- as.data.frame(CC) %>%
    mutate(Category = "CC")
  
  MF <- as.data.frame(MF) %>%
    mutate(Category = "MF")
  
  all_go <- bind_rows(BP, CC, MF)
  
  write.table(
    all_go,
    file = file.path(
      output_dir,
      "GO_enrichment",
      paste0("Williams82_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}