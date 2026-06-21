###############################################################################
## Step 2: Integration and clustering using Harmony
###############################################################################

rm(list = ls())
set.seed(123)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(harmony)

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "1_filter/scn_combine_filter.rds"
output_dir <- "2_integration"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load filtered Seurat object
###############################################################################

scn_combine_filter <- readRDS(input_rds)

DefaultAssay(scn_combine_filter) <- "RNA"

###############################################################################
## 4. Define colors for clusters
###############################################################################

cols <- c(
  "#009933", # 0 cortex
  "#fb90b7", # 1 epidermis
  "#f26500", # 2 xylem
  "#FF007F", # 3 epidermis
  "#A0522D", # 4 meristem
  "#215b63", # 5 pericycle
  "#EE7785", # 6 epidermis
  "#964B00", # 7 meristem
  "#9370DB", # 8 endodermis
  "#6a9c78", # 9 pericycle
  "#c2cfd8", # 10 root hair
  "#9B8281", # 11 xylem
  "#8C3E1C", # 12 meristem
  "#8ee4e4", # 13 phloem
  "#4A672E", # 14 pericycle
  "#00A86B", # 15 cortex
  "#FF7F50", # 16 endodermis
  "#3d6cb9", # 17 syncytium
  "#7a57d1", # 18 endodermis
  "#dddddd", # 19 root hair
  "#FF9900", # 20 xylem
  "#ff4273", # 21 procambium
  "#FFA07A", # 22 scn
  "#1fffff", # 23 phloem
  "#556B2F", # 24 pericycle
  "#00c9b1", # 25 cortex
  "#ff9898", # 26 epidermis
  "#FFB6C1", # 27 scn
  "#a3a1a1"  # 28 root hair
)

names(cols) <- as.character(0:28)

cluster_levels <- as.character(0:28)
cluster_cols <- cols[cluster_levels]

###############################################################################
## 5. Harmony integration
###############################################################################

scn_combine_filter_harmony <- IntegrateLayers(
  object = scn_combine_filter,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)

###############################################################################
## 6. Add sample, library and cultivar metadata
###############################################################################

scn_combine_filter_harmony$library <- scn_combine_filter_harmony$orig.ident

scn_combine_filter_harmony$sample <- case_when(
  grepl("^F1_1$|^F1_2$", scn_combine_filter_harmony$orig.ident) ~ "F1",
  grepl("^F3_1$|^F3_2$", scn_combine_filter_harmony$orig.ident) ~ "F3",
  grepl("^F7_1$|^F7_2$", scn_combine_filter_harmony$orig.ident) ~ "F7",
  grepl("^P1_1$|^P1_2$", scn_combine_filter_harmony$orig.ident) ~ "P1",
  grepl("^P3_1$|^P3_2$", scn_combine_filter_harmony$orig.ident) ~ "P3",
  grepl("^P7_1$|^P7_2$", scn_combine_filter_harmony$orig.ident) ~ "P7",
  grepl("^W1_1$|^W1_2$", scn_combine_filter_harmony$orig.ident) ~ "W1",
  grepl("^W3_1$|^W3_2$", scn_combine_filter_harmony$orig.ident) ~ "W3",
  grepl("^W7_1$|^W7_2$", scn_combine_filter_harmony$orig.ident) ~ "W7",
  TRUE ~ NA_character_
)

scn_combine_filter_harmony$breed <- case_when(
  grepl("^F1_1$|^F1_2$|^F3_1$|^F3_2$|^F7_1$|^F7_2$", scn_combine_filter_harmony$orig.ident) ~ "Forrest",
  grepl("^P1_1$|^P1_2$|^P3_1$|^P3_2$|^P7_1$|^P7_2$", scn_combine_filter_harmony$orig.ident) ~ "PI88788",
  grepl("^W1_1$|^W1_2$|^W3_1$|^W3_2$|^W7_1$|^W7_2$", scn_combine_filter_harmony$orig.ident) ~ "Williams82",
  TRUE ~ NA_character_
)

cell_number_by_sample <- as.data.frame(
  table(scn_combine_filter_harmony$sample)
)

colnames(cell_number_by_sample) <- c("sample", "n_cells")

write.csv(
  cell_number_by_sample,
  file = file.path(output_dir, "cell_number_by_sample.csv"),
  row.names = FALSE
)

cell_number_by_library <- as.data.frame(
  table(scn_combine_filter_harmony$library)
)

colnames(cell_number_by_library) <- c("library", "n_cells")

write.csv(
  cell_number_by_library,
  file = file.path(output_dir, "cell_number_by_library.csv"),
  row.names = FALSE
)

cell_number_by_breed <- as.data.frame(
  table(scn_combine_filter_harmony$breed)
)

colnames(cell_number_by_breed) <- c("breed", "n_cells")

write.csv(
  cell_number_by_breed,
  file = file.path(output_dir, "cell_number_by_breed.csv"),
  row.names = FALSE
)

###############################################################################
## 7. Clustering and UMAP
###############################################################################

resolution_use <- 1
dims_use <- 1:30

scn_combine_filter_harmony <- FindNeighbors(
  object = scn_combine_filter_harmony,
  reduction = "harmony",
  dims = dims_use
)

scn_combine_filter_harmony <- FindClusters(
  object = scn_combine_filter_harmony,
  resolution = resolution_use,
  random.seed = 123,
  verbose = TRUE
)

scn_combine_filter_harmony <- RunUMAP(
  object = scn_combine_filter_harmony,
  reduction = "harmony",
  dims = dims_use,
  seed.use = 123,
  verbose = TRUE
)

scn_combine_filter_harmony$seurat_clusters <- factor(
  scn_combine_filter_harmony$seurat_clusters,
  levels = cluster_levels
)

Idents(scn_combine_filter_harmony) <- "seurat_clusters"

###############################################################################
## 8. Export cluster distribution tables
###############################################################################

cluster_breed_table <- table(
  scn_combine_filter_harmony$seurat_clusters,
  scn_combine_filter_harmony$breed
)

cluster_library_table <- table(
  scn_combine_filter_harmony$seurat_clusters,
  scn_combine_filter_harmony$library
)

cluster_breed_percent <- prop.table(
  cluster_breed_table,
  margin = 2
) * 100

cluster_library_percent <- prop.table(
  cluster_library_table,
  margin = 2
) * 100

write.csv(
  as.data.frame.matrix(cluster_breed_table),
  file = file.path(
    output_dir,
    paste0("cluster_breed_count_resolution_", resolution_use, ".csv")
  )
)

write.csv(
  as.data.frame.matrix(cluster_breed_percent),
  file = file.path(
    output_dir,
    paste0("cluster_breed_percent_resolution_", resolution_use, ".csv")
  )
)

write.csv(
  as.data.frame.matrix(cluster_library_table),
  file = file.path(
    output_dir,
    paste0("cluster_library_count_resolution_", resolution_use, ".csv")
  )
)

write.csv(
  as.data.frame.matrix(cluster_library_percent),
  file = file.path(
    output_dir,
    paste0("cluster_library_percent_resolution_", resolution_use, ".csv")
  )
)

###############################################################################
## 9. Stacked bar plot: cluster distribution by cultivar
###############################################################################

cluster_breed_percent_df <- as.data.frame.matrix(
  cluster_breed_percent
)

cluster_breed_percent_df$Cluster <- rownames(
  cluster_breed_percent_df
)

cluster_breed_long <- cluster_breed_percent_df %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Cultivar",
    values_to = "Percent"
  ) %>%
  mutate(
    Cultivar = recode(
      Cultivar,
      "Williams82" = "Williams 82",
      "PI88788" = "PI 88788",
      "Forrest" = "Forrest"
    ),
    Cultivar = factor(
      Cultivar,
      levels = c("Forrest", "PI 88788", "Williams 82")
    ),
    Cluster = factor(
      Cluster,
      levels = cluster_levels
    )
  )

p_cluster_breed <- ggplot(
  cluster_breed_long,
  aes(
    x = Cultivar,
    y = Percent,
    fill = Cluster
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  scale_fill_manual(
    values = cluster_cols,
    breaks = cluster_levels
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Percent (%)",
    fill = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0("cluster_breed_stacked_bar_resolution_", resolution_use, ".pdf")
  ),
  plot = p_cluster_breed,
  width = 12,
  height = 1.3
)

###############################################################################
## 10. Stacked bar plot: cluster distribution by 18 libraries
###############################################################################

library_order <- c(
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

cluster_library_percent_df <- as.data.frame.matrix(
  cluster_library_percent
)

cluster_library_percent_df$Cluster <- rownames(
  cluster_library_percent_df
)

cluster_library_long <- cluster_library_percent_df %>%
  pivot_longer(
    cols = -Cluster,
    names_to = "Library",
    values_to = "Percent"
  ) %>%
  mutate(
    Library = factor(
      Library,
      levels = rev(library_order)
    ),
    Cluster = factor(
      Cluster,
      levels = cluster_levels
    )
  )

p_cluster_library <- ggplot(
  cluster_library_long,
  aes(
    x = Library,
    y = Percent,
    fill = Cluster
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.7
  ) +
  scale_fill_manual(
    values = cluster_cols,
    breaks = cluster_levels
  ) +
  coord_flip() +
  labs(
    x = "",
    y = "Percent (%)",
    fill = "Cluster"
  ) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 10)
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0("cluster_library_stacked_bar_resolution_", resolution_use, ".pdf")
  ),
  plot = p_cluster_library,
  width = 12,
  height = 5
)

###############################################################################
## 11. UMAP visualization
###############################################################################

p_umap_breed_split <- DimPlot(
  object = scn_combine_filter_harmony,
  reduction = "umap",
  split.by = "breed",
  pt.size = 0.000001,
  raster = FALSE
) +
  scale_color_manual(
    values = cluster_cols
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0("umap_breed_resolution_", resolution_use, ".pdf")
  ),
  plot = p_umap_breed_split,
  width = 18.5,
  height = 5
)

p_umap_cluster <- DimPlot(
  object = scn_combine_filter_harmony,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  shuffle = TRUE,
  label.size = 3,
  alpha = 0.8,
  seed = 10
) +
  scale_color_manual(
    values = cluster_cols
  )

ggsave(
  filename = file.path(
    output_dir,
    paste0("umap_cluster_resolution_", resolution_use, ".pdf")
  ),
  plot = p_umap_cluster,
  width = 8,
  height = 5
)

p_umap_sample <- DimPlot(
  object = scn_combine_filter_harmony,
  reduction = "umap",
  group.by = "sample",
  pt.size = 0.000001,
  raster = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    paste0("umap_sample_resolution_", resolution_use, ".pdf")
  ),
  plot = p_umap_sample,
  width = 8,
  height = 5
)

p_umap_library <- DimPlot(
  object = scn_combine_filter_harmony,
  reduction = "umap",
  group.by = "library",
  pt.size = 0.000001,
  raster = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    paste0("umap_library_resolution_", resolution_use, ".pdf")
  ),
  plot = p_umap_library,
  width = 8,
  height = 5
)

###############################################################################
## 12. Find marker genes for each cluster
###############################################################################

DefaultAssay(scn_combine_filter_harmony) <- "RNA"

scn_combine_filter_harmony <- JoinLayers(
  scn_combine_filter_harmony
)

Idents(scn_combine_filter_harmony) <- "seurat_clusters"

scn_combine_filter_harmony.markers <- FindAllMarkers(
  object = scn_combine_filter_harmony,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.5
)

write.table(
  scn_combine_filter_harmony.markers,
  file = file.path(
    output_dir,
    paste0("scn_all_res", resolution_use, ".markers_fc0.5.txt")
  ),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

###############################################################################
## 13. Save integrated object
###############################################################################

saveRDS(
  scn_combine_filter_harmony,
  file = file.path(
    output_dir,
    paste0("scn_combine_filter_harmony_res", resolution_use, ".rds")
  )
)