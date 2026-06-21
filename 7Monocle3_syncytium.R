###############################################################################
## Step 7: Reconstruction of syncytium developmental trajectories
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
library(ClusterGVis)
library(clusterProfiler)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "3_annotation/scn_combine_filter_harmony_celltype.rds"

output_dir <- "7_syncytium_development"

seurat_dir <- file.path(
  output_dir,
  "seurat_reclustering"
)

monocle_dir <- file.path(
  output_dir,
  "monocle3"
)

Forrest_dir <- file.path(
  monocle_dir,
  "Forrest"
)

PI88788_dir <- file.path(
  monocle_dir,
  "PI88788"
)

Williams82_dir <- file.path(
  monocle_dir,
  "Williams82"
)

branch1_dir <- file.path(
  Forrest_dir,
  "Forrest_branch1"
)

branch2_dir <- file.path(
  Forrest_dir,
  "Forrest_branch2"
)

go_dir <- file.path(
  output_dir,
  "GO_enrichment"
)

go_annotation_file <- "annotation/GOannotation_v4.tsv"
go_info_file <- "annotation/go.tb"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(seurat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(monocle_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(Forrest_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(PI88788_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(Williams82_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(branch1_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(branch2_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load annotated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"

scn_combine_filter_harmony <- JoinLayers(
  scn_combine_filter_harmony
)

###############################################################################
## 4. Extract syncytium cells
###############################################################################

syn <- subset(
  scn_combine_filter_harmony,
  subset = cell_type == "Syncytium"
)

syn <- JoinLayers(syn)

###############################################################################
## 5. Add infection time and batch metadata
###############################################################################

syn$time <- case_when(
  grepl("F1_1|F1_2|P1_1|P1_2|W1_1|W1_2", syn$orig.ident) ~ "1 dpi",
  grepl("F3_1|F3_2|P3_1|P3_2|W3_1|W3_2", syn$orig.ident) ~ "3 dpi",
  grepl("F7_1|F7_2|P7_1|P7_2|W7_1|W7_2", syn$orig.ident) ~ "7 dpi",
  TRUE ~ NA_character_
)

syn$batch <- case_when(
  grepl("F1_1|F3_1|F7_1|P1_1|P3_1|P7_1|W1_1|W3_1|W7_1", syn$orig.ident) ~ "batch1",
  grepl("F1_2|F3_2|F7_2|P1_2|P3_2|P7_2|W1_2|W3_2|W7_2", syn$orig.ident) ~ "batch2",
  TRUE ~ NA_character_
)

###############################################################################
## 6. Re-cluster syncytium cells
###############################################################################

syn <- NormalizeData(syn)

syn <- FindVariableFeatures(
  syn,
  selection.method = "vst",
  nfeatures = 2000
)

syn <- ScaleData(syn)

syn <- RunPCA(
  syn,
  npcs = 50
)

pdf(
  file = file.path(seurat_dir, "syncytium_ElbowPlot.pdf"),
  width = 8,
  height = 6
)

print(
  ElbowPlot(
    syn,
    ndims = 50
  )
)

dev.off()

syn <- FindNeighbors(
  syn,
  dims = 1:20
)

syn <- FindClusters(
  object = syn,
  resolution = 0.1,
  verbose = TRUE
)

syn <- RunUMAP(
  object = syn,
  dims = 1:20,
  verbose = TRUE
)

p <- DimPlot(
  object = syn,
  reduction = "umap",
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 0.8,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_cluster_resolution_0.1.pdf"),
  plot = p,
  width = 6,
  height = 5
)

write.csv(
  as.data.frame(
    table(
      syn$RNA_snn_res.0.1,
      syn$breed
    )
  ),
  file = file.path(seurat_dir, "syncytium_cluster_by_breed.csv"),
  row.names = FALSE
)

write.csv(
  as.data.frame(
    table(
      syn$RNA_snn_res.0.1,
      syn$time
    )
  ),
  file = file.path(seurat_dir, "syncytium_cluster_by_time.csv"),
  row.names = FALSE
)

p <- DimPlot(
  object = syn,
  reduction = "umap",
  group.by = "time",
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 0.8,
  seed = 10
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_cluster_resolution_0.1_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

p <- DimPlot(
  object = syn,
  reduction = "umap",
  split.by = "time",
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 1,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_cluster_resolution_0.1_split_time.pdf"),
  plot = p,
  width = 18,
  height = 5
)

p <- DimPlot(
  object = syn,
  reduction = "umap",
  split.by = "breed",
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 1,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_cluster_resolution_0.1_split_breed.pdf"),
  plot = p,
  width = 18,
  height = 5
)

syn <- JoinLayers(
  syn
)

Idents(syn) <- "RNA_snn_res.0.1"

syn.markers <- FindAllMarkers(
  syn,
  only.pos = TRUE,
  min.pct = 0.25
)

write.table(
  syn.markers,
  file = file.path(seurat_dir, "syncytium_res0.1.markers_fc0.1.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

saveRDS(
  syn,
  file = file.path(seurat_dir, "syn_sub_reclustered.rds")
)

###############################################################################
## 7. Split syncytium cells by cultivar
###############################################################################

Forrest <- subset(
  x = syn,
  subset = breed == "Forrest"
)

PI88788 <- subset(
  x = syn,
  subset = breed == "PI88788"
)

W82 <- subset(
  x = syn,
  subset = breed == "Williams82"
)

Forrest <- JoinLayers(Forrest)
PI88788 <- JoinLayers(PI88788)
W82 <- JoinLayers(W82)

p <- DimPlot(
  object = Forrest,
  reduction = "umap",
  split.by = "time",
  ncol = 1,
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 1,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_Forrest_cluster_resolution_0.1_split_time.pdf"),
  plot = p,
  width = 6,
  height = 10
)

write.csv(
  as.data.frame(
    table(
      Forrest$RNA_snn_res.0.1,
      Forrest$time
    )
  ),
  file = file.path(seurat_dir, "Forrest_cluster_by_time.csv"),
  row.names = FALSE
)

p <- DimPlot(
  object = PI88788,
  reduction = "umap",
  split.by = "time",
  ncol = 1,
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 1,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_PI88788_cluster_resolution_0.1_split_time.pdf"),
  plot = p,
  width = 6,
  height = 10
)

p <- DimPlot(
  object = W82,
  reduction = "umap",
  split.by = "time",
  ncol = 1,
  pt.size = 1,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 1,
  seed = 10
) +
  scale_color_manual(
    values = c("#e1135a", "#62c1ba", "#9489fa")
  )

ggsave(
  filename = file.path(seurat_dir, "syncytium_W82_cluster_resolution_0.1_split_time.pdf"),
  plot = p,
  width = 6,
  height = 10
)

###############################################################################
## 8. Forrest syncytium developmental trajectory
###############################################################################

expr_matrix <- as.matrix(
  LayerData(
    Forrest,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grep("ann1", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- Forrest@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(expr_matrix)
)

rownames(gene_annotation) <- rownames(expr_matrix)

cds_Forrest <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds_Forrest <- preprocess_cds(
  cds_Forrest,
  num_dim = 20
)

cds_Forrest <- align_cds(
  cds_Forrest,
  alignment_group = "batch"
)

cds_Forrest <- reduce_dimension(
  cds_Forrest,
  preprocess_method = "PCA"
)

cds_Forrest <- cluster_cells(
  cds_Forrest,
  reduction_method = "UMAP"
)

p <- plot_cells(
  cds_Forrest,
  reduction_method = "UMAP"
)

ggsave(
  filename = file.path(Forrest_dir, "cds_Forrest_recluster.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_Forrest_embed <- cds_Forrest@int_colData$reducedDims$UMAP

int_Forrest_embed <- Embeddings(
  Forrest,
  reduction = "umap"
)

int_Forrest_embed <- int_Forrest_embed[
  rownames(cds_Forrest_embed),
  ,
  drop = FALSE
]

cds_Forrest@int_colData$reducedDims$UMAP <- int_Forrest_embed

p <- plot_cells(
  cds_Forrest,
  reduction_method = "UMAP",
  color_cells_by = "time"
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(Forrest_dir, "cds_Forrest_umap_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_Forrest <- learn_graph(
  cds_Forrest,
  use_partition = FALSE
)

p <- plot_cells(
  cds_Forrest,
  color_cells_by = "time",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(Forrest_dir, "cds_Forrest_cell_type.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_Forrest <- order_cells(
  cds_Forrest
)

p <- plot_cells(
  cds_Forrest,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "#606470",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  filename = file.path(Forrest_dir, "cds_Forrest_pseudotime.pdf"),
  plot = p,
  width = 7,
  height = 5
)

pr_test_res_Forrest <- graph_test(
  cds_Forrest,
  neighbor_graph = "principal_graph",
  cores = 4
)

pr_deg_ids_Forrest <- rownames(
  subset(
    pr_test_res_Forrest,
    q_value <= 0.01
  )
)

write.table(
  pr_test_res_Forrest,
  file = file.path(Forrest_dir, "Forrest_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  pr_deg_ids_Forrest,
  file = file.path(Forrest_dir, "Forrest_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

pre_pseudotime_matrix <- getFromNamespace(
  "pre_pseudotime_matrix",
  "ClusterGVis"
)

mat_Forrest <- pre_pseudotime_matrix(
  cds_obj = cds_Forrest,
  gene_list = pr_deg_ids_Forrest
)

ck_Forrest <- clusterData(
  obj = mat_Forrest,
  clusterMethod = "kmeans",
  clusterNum = 4
)

Forrest_mark_genes <- character(0)

Forrest_tf_file <- file.path(
  Forrest_dir,
  "Forrest_syncytium_TF.txt"
)

if (file.exists(Forrest_tf_file)) {
  For_tf <- read.table(
    Forrest_tf_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  Forrest_mark_genes <- as.vector(
    For_tf$V1
  )
}

pdf(
  file = file.path(Forrest_dir, "Forrest_pseudotime_pheatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck_Forrest,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = Forrest_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(1, 3, 2)
)

dev.off()

Forrest_heatmap_genes <- as.data.frame(
  ck_Forrest$long.res
)

Forrest_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  Forrest_heatmap_genes$gene
)

write.table(
  Forrest_heatmap_genes,
  file = file.path(Forrest_dir, "Forrest_pseudotime_pheatmap_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 9. Forrest branch1: interactive branch1 selection
###############################################################################

cds_branch1 <- choose_graph_segments(
  cds_Forrest
)

branch1_cells <- colnames(
  cds_branch1
)

cds_branch1 <- cds_Forrest[
  ,
  branch1_cells
]

cds_branch1 <- preprocess_cds(
  cds_branch1,
  num_dim = 20
)

cds_branch1 <- align_cds(
  cds_branch1,
  alignment_group = "batch"
)

cds_branch1 <- reduce_dimension(
  cds_branch1,
  preprocess_method = "PCA",
  reduction_method = "UMAP"
)

branch1_embed <- Embeddings(
  Forrest,
  reduction = "umap"
)

branch1_embed <- branch1_embed[
  colnames(cds_branch1),
  ,
  drop = FALSE
]

reducedDims(cds_branch1)$UMAP <- branch1_embed

cds_branch1 <- cluster_cells(
  cds_branch1,
  reduction_method = "UMAP"
)

cds_branch1 <- learn_graph(
  cds_branch1,
  use_partition = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 15
  )
)

cds_branch1 <- order_cells(
  cds_branch1
)

p <- plot_cells(
  cds_branch1,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "#606470",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  filename = file.path(branch1_dir, "cds_Forrest_branch1_pseudotime.pdf"),
  plot = p,
  width = 7,
  height = 5
)

p <- plot_cells(
  cds_branch1,
  color_cells_by = "time",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(branch1_dir, "cds_Forrest_branch1_dpi.pdf"),
  plot = p,
  width = 6,
  height = 5
)

pr_test_res_Forrest_branch1 <- graph_test(
  cds_branch1,
  neighbor_graph = "principal_graph",
  cores = 4
)

pr_deg_ids_Forrest_branch1 <- rownames(
  subset(
    pr_test_res_Forrest_branch1,
    q_value <= 0.01
  )
)

write.table(
  pr_test_res_Forrest_branch1,
  file = file.path(branch1_dir, "Forrest_branch1_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

pre_pseudotime_matrix_cds <- getFromNamespace(
  "pre_pseudotime_matrix",
  "ClusterGVis"
)

mat_Forrest_branch1 <- pre_pseudotime_matrix_cds(
  cds_obj = cds_branch1,
  gene_list = pr_deg_ids_Forrest_branch1
)

ck_Forrest_branch1 <- clusterData(
  obj = mat_Forrest_branch1,
  clusterMethod = "kmeans",
  clusterNum = 3
)

Forrest_branch1_genes <- ck_Forrest_branch1$long.res %>%
  select(cluster, gene) %>%
  distinct() %>%
  arrange(cluster, gene)

Forrest_branch1_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  Forrest_branch1_genes$gene
)

write.table(
  Forrest_branch1_genes,
  file = file.path(branch1_dir, "Forrest_branch1_pseudotime_pheatmap_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

Forrest_branch1_mark_genes <- character(0)

Forrest_branch1_tf_file <- file.path(
  branch1_dir,
  "Forrest_branch1_TF.txt"
)

if (file.exists(Forrest_branch1_tf_file)) {
  For_tf_branch1 <- read.table(
    Forrest_branch1_tf_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  Forrest_branch1_mark_genes <- as.vector(
    For_tf_branch1$V1
  )
}

pdf(
  file = file.path(branch1_dir, "Forrest_branch1_pseudotime_pheatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck_Forrest_branch1,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = Forrest_branch1_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(3, 1, 2)
)

dev.off()

###############################################################################
## 10. Forrest branch2: interactive branch2 selection
###############################################################################

cds_branch2 <- choose_graph_segments(
  cds_Forrest
)

branch2_cells <- colnames(
  cds_branch2
)

cds_branch2 <- cds_Forrest[
  ,
  branch2_cells
]

cds_branch2 <- preprocess_cds(
  cds_branch2,
  num_dim = 20
)

cds_branch2 <- align_cds(
  cds_branch2,
  alignment_group = "batch"
)

cds_branch2 <- reduce_dimension(
  cds_branch2,
  preprocess_method = "PCA",
  reduction_method = "UMAP"
)

branch2_embed <- Embeddings(
  Forrest,
  reduction = "umap"
)

branch2_embed <- branch2_embed[
  colnames(cds_branch2),
  ,
  drop = FALSE
]

reducedDims(cds_branch2)$UMAP <- branch2_embed

cds_branch2 <- cluster_cells(
  cds_branch2,
  reduction_method = "UMAP"
)

cds_branch2 <- learn_graph(
  cds_branch2,
  use_partition = FALSE,
  learn_graph_control = list(
    minimal_branch_len = 15
  )
)

cds_branch2 <- order_cells(
  cds_branch2
)

p <- plot_cells(
  cds_branch2,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "#606470",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  filename = file.path(branch2_dir, "cds_Forrest_branch2_pseudotime.pdf"),
  plot = p,
  width = 7,
  height = 5
)

p <- plot_cells(
  cds_branch2,
  color_cells_by = "time",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(branch2_dir, "cds_Forrest_branch2_dpi.pdf"),
  plot = p,
  width = 6,
  height = 5
)

pr_test_res_Forrest_branch2 <- graph_test(
  cds_branch2,
  neighbor_graph = "principal_graph",
  cores = 4
)

pr_deg_ids_Forrest_branch2 <- rownames(
  subset(
    pr_test_res_Forrest_branch2,
    q_value <= 0.01
  )
)

write.table(
  pr_test_res_Forrest_branch2,
  file = file.path(branch2_dir, "Forrest_branch2_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

pre_pseudotime_matrix_cds <- getFromNamespace(
  "pre_pseudotime_matrix",
  "ClusterGVis"
)

mat_Forrest_branch2 <- pre_pseudotime_matrix_cds(
  cds_obj = cds_branch2,
  gene_list = pr_deg_ids_Forrest_branch2
)

ck_Forrest_branch2 <- clusterData(
  obj = mat_Forrest_branch2,
  clusterMethod = "kmeans",
  clusterNum = 3
)

Forrest_branch2_genes <- ck_Forrest_branch2$long.res %>%
  select(cluster, gene) %>%
  distinct() %>%
  arrange(cluster, gene)

Forrest_branch2_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  Forrest_branch2_genes$gene
)

write.table(
  Forrest_branch2_genes,
  file = file.path(branch2_dir, "Forrest_branch2_pseudotime_pheatmap_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

Forrest_branch2_mark_genes <- character(0)

Forrest_branch2_tf_file <- file.path(
  branch2_dir,
  "Forrest_branch2_TF.txt"
)

if (file.exists(Forrest_branch2_tf_file)) {
  For_tf_branch2 <- read.table(
    Forrest_branch2_tf_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  Forrest_branch2_mark_genes <- as.vector(
    For_tf_branch2$V1
  )
}

pdf(
  file = file.path(branch2_dir, "Forrest_branch2_pseudotime_pheatmap.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck_Forrest_branch2,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = Forrest_branch2_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(2, 3, 1)
)

dev.off()

###############################################################################
## 11. Time distribution of Forrest branch1 and branch2
###############################################################################

branch1_df <- data.frame(
  cell = colnames(cds_branch1),
  branch = "Branch1",
  time = colData(cds_branch1)$time
)

branch2_df <- data.frame(
  cell = colnames(cds_branch2),
  branch = "Branch2",
  time = colData(cds_branch2)$time
)

branch_time_df <- bind_rows(
  branch1_df,
  branch2_df
)

branch_time_stat <- branch_time_df %>%
  group_by(branch, time) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(branch) %>%
  mutate(percent = 100 * count / sum(count)) %>%
  ungroup()

write.table(
  as.data.frame(branch_time_stat),
  file = file.path(Forrest_dir, "Forrest_two_branch_dpi_origin.txt"),
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)

p <- ggplot(
  branch_time_stat,
  aes(x = branch, y = percent, fill = time)
) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  ) +
  labs(
    x = "",
    y = "Proportion (%)",
    fill = "Time"
  ) +
  theme_classic()

ggsave(
  filename = file.path(Forrest_dir, "Forrest_two_branch_dpi_origin.pdf"),
  plot = p,
  width = 4,
  height = 7
)

###############################################################################
## 12. PI88788 syncytium developmental trajectory
###############################################################################

expr_matrix <- as.matrix(
  LayerData(
    PI88788,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grep("ann1", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- PI88788@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(expr_matrix)
)

rownames(gene_annotation) <- rownames(expr_matrix)

cds_PI88788 <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds_PI88788 <- preprocess_cds(
  cds_PI88788,
  num_dim = 20
)

cds_PI88788 <- align_cds(
  cds_PI88788,
  alignment_group = "time"
)

cds_PI88788 <- reduce_dimension(
  cds_PI88788,
  preprocess_method = "PCA"
)

cds_PI88788 <- cluster_cells(
  cds_PI88788,
  reduction_method = "UMAP"
)

p <- plot_cells(
  cds_PI88788,
  reduction_method = "UMAP"
)

ggsave(
  filename = file.path(PI88788_dir, "cds_PI88788_recluster.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_PI88788_embed <- cds_PI88788@int_colData$reducedDims$UMAP

int_PI88788_embed <- Embeddings(
  PI88788,
  reduction = "umap"
)

int_PI88788_embed <- int_PI88788_embed[
  rownames(cds_PI88788_embed),
  ,
  drop = FALSE
]

cds_PI88788@int_colData$reducedDims$UMAP <- int_PI88788_embed

p <- plot_cells(
  cds_PI88788,
  reduction_method = "UMAP",
  color_cells_by = "time"
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(PI88788_dir, "cds_PI88788_umap_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_PI88788 <- learn_graph(
  cds_PI88788,
  use_partition = FALSE
)

p <- plot_cells(
  cds_PI88788,
  color_cells_by = "time",
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(PI88788_dir, "cds_PI88788_cell_type.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_PI88788 <- order_cells(
  cds_PI88788
)

p <- plot_cells(
  cds_PI88788,
  color_cells_by = "pseudotime",
  trajectory_graph_color = "#606470",
  show_trajectory_graph = TRUE,
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  filename = file.path(PI88788_dir, "cds_PI88788_pseudotime.pdf"),
  plot = p,
  width = 7,
  height = 5
)

p <- plot_cells(
  cds_PI88788,
  color_cells_by = "time",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(PI88788_dir, "cds_PI88788_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

pr_test_res_PI88788 <- graph_test(
  cds_PI88788,
  neighbor_graph = "principal_graph",
  cores = 4
)

pr_deg_ids_PI88788 <- rownames(
  subset(
    pr_test_res_PI88788,
    q_value <= 0.01
  )
)

write.table(
  pr_test_res_PI88788,
  file = file.path(PI88788_dir, "PI88788_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  pr_deg_ids_PI88788,
  file = file.path(PI88788_dir, "PI88788_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

mat_PI88788 <- pre_pseudotime_matrix(
  cds_obj = cds_PI88788,
  gene_list = pr_deg_ids_PI88788
)

ck_PI88788 <- clusterData(
  obj = mat_PI88788,
  clusterMethod = "kmeans",
  clusterNum = 3
)

PI88788_mark_genes <- character(0)

PI88788_tf_file <- file.path(
  PI88788_dir,
  "PI88788_syncytium_TF.txt"
)

if (file.exists(PI88788_tf_file)) {
  PI_tf <- read.table(
    PI88788_tf_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  PI88788_mark_genes <- as.vector(
    PI_tf$V1
  )
}

pdf(
  file = file.path(PI88788_dir, "PI88788_pseudotime_pheatmap_1.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck_PI88788,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = PI88788_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(3, 1, 2)
)

dev.off()

PI88788_heatmap_genes <- as.data.frame(
  ck_PI88788$long.res
)

PI88788_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  PI88788_heatmap_genes$gene
)

write.table(
  PI88788_heatmap_genes,
  file = file.path(PI88788_dir, "PI88788_pseudotime_pheatmap_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 13. Williams82 syncytium developmental trajectory
###############################################################################

expr_matrix <- as.matrix(
  LayerData(
    W82,
    assay = "RNA",
    layer = "counts"
  )
)

expr_matrix <- as.matrix(
  expr_matrix[
    grep("ann1", rownames(expr_matrix)),
    ,
    drop = FALSE
  ]
)

cell_metadata <- W82@meta.data

gene_annotation <- data.frame(
  gene_short_name = rownames(expr_matrix)
)

rownames(gene_annotation) <- rownames(expr_matrix)

cds_W82 <- new_cell_data_set(
  expr_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_annotation
)

cds_W82 <- preprocess_cds(
  cds_W82,
  num_dim = 20
)

cds_W82 <- align_cds(
  cds_W82,
  alignment_group = "time"
)

cds_W82 <- reduce_dimension(
  cds_W82,
  preprocess_method = "PCA"
)

cds_W82 <- cluster_cells(
  cds_W82,
  reduction_method = "UMAP"
)

p <- plot_cells(
  cds_W82,
  reduction_method = "UMAP"
)

ggsave(
  filename = file.path(Williams82_dir, "cds_W82_recluster.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_W82_embed <- cds_W82@int_colData$reducedDims$UMAP

int_W82_embed <- Embeddings(
  W82,
  reduction = "umap"
)

int_W82_embed <- int_W82_embed[
  rownames(cds_W82_embed),
  ,
  drop = FALSE
]

cds_W82@int_colData$reducedDims$UMAP <- int_W82_embed

p <- plot_cells(
  cds_W82,
  reduction_method = "UMAP",
  color_cells_by = "time"
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(Williams82_dir, "cds_W82_umap_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_W82 <- learn_graph(
  cds_W82,
  use_partition = FALSE
)

p <- plot_cells(
  cds_W82,
  color_cells_by = "time",
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(Williams82_dir, "cds_W82_cell_type.pdf"),
  plot = p,
  width = 6,
  height = 5
)

cds_W82 <- order_cells(
  cds_W82
)

p <- plot_cells(
  cds_W82,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_cell_groups = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE,
  graph_label_size = 1.5,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
)

ggsave(
  filename = file.path(Williams82_dir, "cds_W82_pseudotime.pdf"),
  plot = p,
  width = 7,
  height = 5
)

p <- plot_cells(
  cds_W82,
  color_cells_by = "time",
  show_trajectory_graph = TRUE,
  trajectory_graph_color = "#606470",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  group_label_size = 4,
  cell_size = 1,
  min_expr = 5
) +
  scale_color_manual(
    values = c(
      "1 dpi" = "#FFCC00",
      "3 dpi" = "#009999",
      "7 dpi" = "#CC3333"
    )
  )

ggsave(
  filename = file.path(Williams82_dir, "cds_W82_time.pdf"),
  plot = p,
  width = 6,
  height = 5
)

pr_test_res_W82 <- graph_test(
  cds_W82,
  neighbor_graph = "principal_graph",
  cores = 4
)

pr_deg_ids_W82 <- rownames(
  subset(
    pr_test_res_W82,
    q_value <= 0.01
  )
)

write.table(
  pr_test_res_W82,
  file = file.path(Williams82_dir, "W82_graph_test_all_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = TRUE
)

write.table(
  pr_deg_ids_W82,
  file = file.path(Williams82_dir, "W82_trajectory_genes.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

mat_W82 <- pre_pseudotime_matrix(
  cds_obj = cds_W82,
  gene_list = pr_deg_ids_W82
)

ck_W82 <- clusterData(
  obj = mat_W82,
  clusterMethod = "kmeans",
  clusterNum = 3
)

W82_mark_genes <- character(0)

W82_tf_file <- file.path(
  Williams82_dir,
  "W82_syncytium_TF.txt"
)

if (file.exists(W82_tf_file)) {
  W82_tf <- read.table(
    W82_tf_file,
    header = FALSE,
    stringsAsFactors = FALSE
  )
  
  W82_mark_genes <- as.vector(
    W82_tf$V1
  )
}

pdf(
  file = file.path(Williams82_dir, "W82_pseudotime_pheatmap_1.pdf"),
  height = 10,
  width = 8,
  onefile = FALSE
)

visCluster(
  object = ck_W82,
  plotType = "heatmap",
  addSampleAnno = FALSE,
  markGenes = W82_mark_genes,
  genesGp = c("italic", 6, NA),
  clusterOrder = c(1, 2, 3)
)

dev.off()

W82_heatmap_genes <- as.data.frame(
  ck_W82$long.res
)

W82_heatmap_genes$gene_clean <- sub(
  "^ann1\\.",
  "",
  W82_heatmap_genes$gene
)

write.table(
  W82_heatmap_genes,
  file = file.path(Williams82_dir, "W82_pseudotime_pheatmap_gene.txt"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

###############################################################################
## 14. Syncytium proportion across time points
###############################################################################

syn_breed_time <- table(
  syn$breed,
  syn$time
)

syn_breed_time_ratio <- round(
  syn_breed_time / rowSums(syn_breed_time) * 100,
  2
)

write.csv(
  as.data.frame.matrix(syn_breed_time),
  file = file.path(output_dir, "syncytium_cell_number_by_breed_time.csv")
)

write.csv(
  as.data.frame.matrix(syn_breed_time_ratio),
  file = file.path(output_dir, "syncytium_cell_ratio_by_breed_time.csv")
)

pdf(
  file = file.path(output_dir, "syn_breed_time.pdf"),
  width = 5.8,
  height = 3.3
)

par(
  mar = c(5, 5, 4, 5),
  xpd = TRUE
)

barplot(
  t(syn_breed_time_ratio),
  col = c("#FFCC00", "#009999", "#CC3333"),
  xlab = "Cell proportion (%)",
  horiz = TRUE,
  las = 1
)

legend(
  "top",
  inset = -0.4,
  legend = c("1 dpi", "3 dpi", "7 dpi"),
  pch = 15,
  col = c("#FFCC00", "#009999", "#CC3333"),
  title = "Time point",
  ncol = 3
)

dev.off()

###############################################################################
## 15. GO enrichment for pseudotime heatmap genes
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
## 15.1 Forrest branch1 GO enrichment
###############################################################################

Forrest_branch1_heatmap_gene_table <- read.table(
  file.path(branch1_dir, "Forrest_branch1_pseudotime_pheatmap_gene.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

for (module_id in sort(unique(Forrest_branch1_heatmap_gene_table$cluster))) {
  
  marker <- unique(
    Forrest_branch1_heatmap_gene_table$gene_clean[
      Forrest_branch1_heatmap_gene_table$cluster == module_id
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
  
  all_go <- bind_rows(
    BP,
    CC,
    MF
  )
  
  write.table(
    all_go,
    file = file.path(
      go_dir,
      paste0("Forrest_branch1_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 15.2 Forrest branch2 GO enrichment
###############################################################################

Forrest_branch2_heatmap_gene_table <- read.table(
  file.path(branch2_dir, "Forrest_branch2_pseudotime_pheatmap_gene.txt"),
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

for (module_id in sort(unique(Forrest_branch2_heatmap_gene_table$cluster))) {
  
  marker <- unique(
    Forrest_branch2_heatmap_gene_table$gene_clean[
      Forrest_branch2_heatmap_gene_table$cluster == module_id
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
  
  all_go <- bind_rows(
    BP,
    CC,
    MF
  )
  
  write.table(
    all_go,
    file = file.path(
      go_dir,
      paste0("Forrest_branch2_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 15.3 PI88788 GO enrichment
###############################################################################

PI88788_heatmap_gene_table <- read.table(
  file.path(PI88788_dir, "PI88788_pseudotime_pheatmap_gene.txt"),
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
  
  all_go <- bind_rows(
    BP,
    CC,
    MF
  )
  
  write.table(
    all_go,
    file = file.path(
      go_dir,
      paste0("PI88788_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 15.4 Williams82 GO enrichment
###############################################################################

Williams82_heatmap_gene_table <- read.table(
  file.path(Williams82_dir, "W82_pseudotime_pheatmap_gene.txt"),
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
  
  all_go <- bind_rows(
    BP,
    CC,
    MF
  )
  
  write.table(
    all_go,
    file = file.path(
      go_dir,
      paste0("Williams82_pseudotime_module_", module_id, "_GO_all.txt")
    ),
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
  )
}

###############################################################################
## 16. Save final object
###############################################################################

saveRDS(
  syn,
  file = file.path(output_dir, "syncytium_reclustered_object.rds")
)