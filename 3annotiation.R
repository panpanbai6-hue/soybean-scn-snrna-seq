###############################################################################
## Step 3: Cell type annotation
###############################################################################

rm(list = ls())
set.seed(123)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(MetBrewer)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "2_integration/scn_combine_filter_harmony_res1.rds"
output_dir <- "3_annotation"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load integrated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"
Idents(scn_combine_filter_harmony) <- "seurat_clusters"

###############################################################################
## 4. Define marker genes for cell type annotation
###############################################################################

gene_all <- c(
  ## Xylem
  "Glyma.04G063800",
  "Glyma.06G065000",
  "Glyma.13G334500",
  "Glyma.15G040000",
  
  ## Syncytium
  "Glyma.18G022500",
  
  ## Cortex
  "Glyma.06G143000",
  "Glyma.15G169100",
  "Glyma.19G124500",
  
  ## Pericycle
  "Glyma.02G003700",
  "Glyma.05G023700",
  "Glyma.11G078300",
  
  ## Procambium
  "Glyma.17G235000",
  "Glyma.17G180400",
  
  ## Root hair
  "Glyma.07G130800",
  "Glyma.02G149100",
  "Glyma.17G133100",
  "Glyma.13G293500",
  
  ## Epidermis
  "Glyma.06G259400",
  "Glyma.09G099900",
  "Glyma.10G070200",
  "Glyma.01G156200",
  
  ## Endodermis
  "Glyma.14G218700",
  "Glyma.16G106800",
  "Glyma.11G008000",
  "Glyma.12G030300",
  "Glyma.12G196600",
  
  ## Meristem
  "Glyma.08G344200",
  "Glyma.11G101500",
  "Glyma.12G064300",
  "Glyma.08G346500",
  "Glyma.08G053600",
  "Glyma.14G086300",
  "Glyma.19G239200",
  
  ## Phloem
  "Glyma.05G216000",
  "Glyma.15G274200",
  "Glyma.11G243100",
  "Glyma.12G154300",
  "Glyma.08G308700",
  "Glyma.07G174600",
  "Glyma.15G065400",
  "Glyma.12G213200",
  "Glyma.19G231300"
)

###############################################################################
## 5. Generate marker gene dotplot for soybean nuclei
###############################################################################

scn_combine_filter_harmony_soybean <- subset(
  scn_combine_filter_harmony,
  subset = species == "soybean"
)

Idents(scn_combine_filter_harmony_soybean) <- "seurat_clusters"

marker_features <- paste0("ann1.", gene_all)

marker_features <- intersect(
  marker_features,
  rownames(scn_combine_filter_harmony_soybean)
)


p_marker_dotplot <- DotPlot(
  object = scn_combine_filter_harmony_soybean,
  features = marker_features,
  cluster.idents = TRUE,
  dot.scale = 6,
  scale.min = 0,
  col.max = 2.7,
  col.min = -2
) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  scale_x_discrete(
    labels = function(x) sub("^ann1\\.", "", x)
  ) +
  scale_color_gradientn(
    colors = MetBrewer::met.brewer("Benedictus", direction = -1)
  )

ggsave(
  filename = file.path(output_dir, "markergene_dotplot_cluster.pdf"),
  plot = p_marker_dotplot,
  width = 11,
  height = 6
)

###############################################################################
## 6. Assign cell types to clusters
###############################################################################

cluster_to_celltype <- c(
  "0"  = "Cortex",
  "1"  = "Epidermis",
  "2"  = "Unknown",
  "3"  = "Epidermis",
  "4"  = "Meristem",
  "5"  = "Pericycle",
  "6"  = "Epidermis",
  "7"  = "Meristem",
  "8"  = "Endodermis",
  "9"  = "Pericycle",
  "10" = "Root hair",
  "11" = "Xylem",
  "12" = "Meristem",
  "13" = "Phloem",
  "14" = "Pericycle",
  "15" = "Cortex",
  "16" = "Endodermis",
  "17" = "Syncytium",
  "18" = "Endodermis",
  "19" = "Unknown",
  "20" = "Xylem",
  "21" = "Procambium",
  "22" = "SCN",
  "23" = "Phloem",
  "24" = "Pericycle",
  "25" = "Cortex",
  "26" = "Epidermis",
  "27" = "SCN",
  "28" = "Root hair"
)

cluster_ids <- as.character(
  scn_combine_filter_harmony$seurat_clusters
)

scn_combine_filter_harmony$cell_type <- cluster_to_celltype[cluster_ids]

cluster_celltype_table <- as.data.frame(
  table(
    scn_combine_filter_harmony$seurat_clusters,
    scn_combine_filter_harmony$cell_type
  )
)

colnames(cluster_celltype_table) <- c(
  "cluster",
  "cell_type",
  "n_cells"
)

write.csv(
  cluster_celltype_table,
  file = file.path(output_dir, "cluster_celltype_table.csv"),
  row.names = FALSE
)

cell_number_by_celltype <- as.data.frame(
  table(scn_combine_filter_harmony$cell_type)
)

colnames(cell_number_by_celltype) <- c(
  "cell_type",
  "n_cells"
)

write.csv(
  cell_number_by_celltype,
  file = file.path(output_dir, "cell_number_by_celltype.csv"),
  row.names = FALSE
)

###############################################################################
## 7. UMAP visualization by annotated cell type
###############################################################################

celltype_cols <- c(
  "Cortex"     = "#1fab89",
  "Endodermis" = "#9370DB",
  "Epidermis"  = "#EE7785",
  "Meristem"   = "#964B00",
  "Pericycle"  = "#556B2F",
  "Phloem"     = "#9ddcdc",
  "Procambium" = "#C03B4C",
  "Root hair"  = "#a3a1a1",
  "SCN"        = "#FFA07A",
  "Syncytium"  = "#3d6cb9",
  "Unknown"    = "#4f5e7f",
  "Xylem"      = "#ec610a"
)

p_celltype_umap <- DimPlot(
  object = scn_combine_filter_harmony,
  reduction = "umap",
  group.by = "cell_type",
  pt.size = 0.00001,
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 0.3,
  seed = 10
) +
  scale_color_manual(values = celltype_cols)

ggsave(
  filename = file.path(output_dir, "umap_celltype_resolution_1.pdf"),
  plot = p_celltype_umap,
  width = 8,
  height = 5
)

###############################################################################
## 8. Rename identities by cell type
###############################################################################

Idents(scn_combine_filter_harmony) <- "seurat_clusters"

new_cluster_ids <- cluster_to_celltype[
  levels(scn_combine_filter_harmony)
]

names(new_cluster_ids) <- levels(
  scn_combine_filter_harmony
)

scn_combine_filter_harmony <- RenameIdents(
  object = scn_combine_filter_harmony,
  new_cluster_ids
)

scn_combine_filter_harmony$cell_type_ident <- Idents(
  scn_combine_filter_harmony
)

###############################################################################
## 9. Identify marker genes for annotated cell types
###############################################################################

DefaultAssay(scn_combine_filter_harmony) <- "RNA"

scn_combine_filter_harmony <- JoinLayers(
  scn_combine_filter_harmony
)

Idents(scn_combine_filter_harmony) <- "cell_type"

scn_combine_filter_harmony.markers <- FindAllMarkers(
  object = scn_combine_filter_harmony,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

write.table(
  scn_combine_filter_harmony.markers,
  file = file.path(output_dir, "scn_all_res1_markers_celltype.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

###############################################################################
## 10. Save annotated Seurat object
###############################################################################

saveRDS(
  scn_combine_filter_harmony,
  file = file.path(output_dir, "scn_combine_filter_harmony_celltype.rds")
)