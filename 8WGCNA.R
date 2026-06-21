###############################################################################
## Step 8: hdWGCNA analysis of F3 syncytium nuclei
###############################################################################

rm(list = ls())
set.seed(12345)

###############################################################################
## 1. Load required R packages
###############################################################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(clusterProfiler)
library(igraph)
library(tidygraph)
library(ggraph)

theme_set(theme_cowplot())

###############################################################################
## 2. Set input and output paths
###############################################################################

input_rds <- "3_annotation/scn_combine_filter_harmony_celltype.rds"

go_annotation_file <- "annotation/GOannotation_v4.tsv"
go_info_file <- "annotation/go.tb"

output_dir <- "8_hdWGCNA_syncytium"
hdwgcna_dir <- file.path(output_dir, "hdWGCNA")
plot_dir <- file.path(output_dir, "plots")
go_dir <- file.path(output_dir, "GO_enrichment")
network_dir <- file.path(output_dir, "gene_networks")
manual_list_dir <- file.path(output_dir, "manual_gene_lists")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(hdwgcna_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(network_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(manual_list_dir, recursive = TRUE, showWarnings = FALSE)

###############################################################################
## 3. Load GO annotation
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
## 4. Load annotated Seurat object
###############################################################################

scn_combine_filter_harmony <- readRDS(input_rds)

DefaultAssay(scn_combine_filter_harmony) <- "RNA"
scn_combine_filter_harmony <- JoinLayers(scn_combine_filter_harmony)

required_meta <- c("cell_type", "sample", "orig.ident")

missing_meta <- setdiff(
  required_meta,
  colnames(scn_combine_filter_harmony@meta.data)
)

###############################################################################
## 5. Extract F3 syncytium nuclei
###############################################################################

F3 <- subset(
  scn_combine_filter_harmony,
  subset = sample == "F3"
)

F3_syn <- subset(
  F3,
  subset = cell_type == "Syncytium"
)

F3_syn <- JoinLayers(F3_syn)

###############################################################################
## 6. Setup hdWGCNA gene set
###############################################################################

F3_syn <- SetupForWGCNA(
  F3_syn,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "soybean"
)

wgcna_genes <- F3_syn@misc$soybean$wgcna_genes

write.table(
  data.frame(gene = wgcna_genes),
  file = file.path(hdwgcna_dir, "wgcna_genes_fraction_0.05.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

message("Number of genes used for hdWGCNA: ", length(wgcna_genes))

###############################################################################
## 7. Construct metacells
###############################################################################

F3_syn <- MetacellsByGroups(
  seurat_obj = F3_syn,
  reduction = "harmony",
  k = 25,
  max_shared = 10,
  min_cells = 70,
  dims = 1:15
)

F3_syn <- NormalizeMetacells(F3_syn)

F3_syn <- ScaleMetacells(
  F3_syn,
  features = VariableFeatures(F3_syn)
)

F3_syn <- RunPCAMetacells(
  F3_syn,
  features = VariableFeatures(F3_syn)
)

F3_syn <- RunUMAPMetacells(
  F3_syn,
  dim = 1:30
)

p_metacell <- DimPlotMetacells(F3_syn) +
  umap_theme() +
  ggtitle("F3 Syncytium metacells")

ggsave(
  filename = file.path(plot_dir, "F3_syncytium_DimPlotMetacells.pdf"),
  plot = p_metacell,
  width = 7.5,
  height = 5
)

###############################################################################
## 8. Prepare expression matrix for hdWGCNA
###############################################################################

F3_syn <- SetDatExpr(
  F3_syn,
  assay = "RNA",
  slot = "data"
)

###############################################################################
## 9. Test soft-thresholding powers
###############################################################################

F3_syn <- TestSoftPowers(
  F3_syn,
  networkType = "signed"
)

plot_list <- PlotSoftPowers(F3_syn)

p_soft_power <- wrap_plots(plot_list, ncol = 2)

ggsave(
  filename = file.path(plot_dir, "soft_power_threshold.pdf"),
  plot = p_soft_power,
  width = 8,
  height = 8
)

power_table <- GetPowerTable(F3_syn)

write.table(
  power_table,
  file = file.path(hdwgcna_dir, "soft_power_table.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

###############################################################################
## 10. Construct co-expression network
###############################################################################

soft_power_use <- 4

F3_syn <- ConstructNetwork(
  F3_syn,
  soft_power = soft_power_use,
  setDatExpr = FALSE,
  tom_name = "Syncytium",
  overwrite_tom = TRUE
)

pdf(
  file = file.path(plot_dir, "Syncytium_hdWGCNA_Dendrogram.pdf"),
  width = 5,
  height = 4
)

PlotDendrogram(
  F3_syn,
  main = "Syncytium hdWGCNA Dendrogram"
)

dev.off()

###############################################################################
## 11. Compute module eigengenes and module connectivity
###############################################################################

F3_syn <- ScaleData(
  F3_syn,
  features = VariableFeatures(F3_syn)
)

F3_syn <- ModuleEigengenes(F3_syn)

hMEs <- GetMEs(F3_syn)
MEs <- GetMEs(F3_syn, harmonized = FALSE)

write.table(
  hMEs,
  file = file.path(hdwgcna_dir, "harmonized_module_eigengenes.txt"),
  quote = FALSE,
  sep = "\t"
)

write.table(
  MEs,
  file = file.path(hdwgcna_dir, "module_eigengenes.txt"),
  quote = FALSE,
  sep = "\t"
)

F3_syn <- ModuleConnectivity(F3_syn)

F3_syn <- ResetModuleNames(
  F3_syn,
  new_name = "Syncytium-M"
)

p_kme <- PlotKMEs(F3_syn, ncol = 5)

ggsave(
  filename = file.path(plot_dir, "module_kME_rank_plot.pdf"),
  plot = p_kme,
  width = 12,
  height = 8
)

###############################################################################
## 12. Export modules and hub genes
###############################################################################

modules <- GetModules(F3_syn) %>%
  subset(module != "grey")

write.table(
  modules,
  file = file.path(hdwgcna_dir, "modules.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

module_count <- as.data.frame(table(modules$module))

write.csv(
  module_count,
  file = file.path(hdwgcna_dir, "module_gene_number.csv"),
  row.names = FALSE
)

hub_df <- GetHubGenes(
  F3_syn,
  n_hubs = 80
)

write.table(
  hub_df,
  file = file.path(hdwgcna_dir, "HubGene_top80.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

###############################################################################
## 13. Extract TOM matrix and Syncytium-M2 module genes
###############################################################################

TOM <- GetTOM(F3_syn)

get_module_genes <- function(module_table, module_name) {
  gene_col_candidates <- c("gene_name", "gene", "gene_id")
  gene_col <- intersect(gene_col_candidates, colnames(module_table))
  
  if (length(gene_col) > 0) {
    genes <- module_table[module_table$module == module_name, gene_col[1]]
  } else {
    genes <- rownames(module_table[module_table$module == module_name, , drop = FALSE])
  }
  
  unique(as.character(genes))
}

Msyn_genes <- get_module_genes(
  module_table = modules,
  module_name = "Syncytium-M2"
)

Msyn_genes <- intersect(Msyn_genes, rownames(TOM))

if (length(Msyn_genes) == 0) {
  stop("No genes were found for module Syncytium-M2.")
}

TOM_syn <- TOM[Msyn_genes, Msyn_genes]

write.table(
  data.frame(gene = sub("^ann1\\.", "", Msyn_genes)),
  file = file.path(hdwgcna_dir, "Syncytium_M2_genes.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

###############################################################################
## 14. GO enrichment for Syncytium-M2 genes
###############################################################################

Msyn_genes_GO <- sub("^ann1\\.", "", Msyn_genes)

BP <- enricher(
  Msyn_genes_GO,
  TERM2GENE = GOannotation[
    GOannotation$level == "biological_process",
    c(2, 1)
  ],
  TERM2NAME = GOinfo[, 1:2],
  pAdjustMethod = "BH"
)

CC <- enricher(
  Msyn_genes_GO,
  TERM2GENE = GOannotation[
    GOannotation$level == "cellular_component",
    c(2, 1)
  ],
  TERM2NAME = GOinfo[, 1:2],
  pAdjustMethod = "BH"
)

MF <- enricher(
  Msyn_genes_GO,
  TERM2GENE = GOannotation[
    GOannotation$level == "molecular_function",
    c(2, 1)
  ],
  TERM2NAME = GOinfo[, 1:2],
  pAdjustMethod = "BH"
)

BP <- as.data.frame(BP)
CC <- as.data.frame(CC)
MF <- as.data.frame(MF)

if (nrow(BP) > 0) BP$Category <- "BP"
if (nrow(CC) > 0) CC$Category <- "CC"
if (nrow(MF) > 0) MF$Category <- "MF"

Msyn_GO <- bind_rows(BP, CC, MF)

write.table(
  Msyn_GO,
  file = file.path(go_dir, "Syncytium_M2_GO_all.txt"),
  row.names = FALSE,
  sep = "\t",
  quote = FALSE
)

###############################################################################
## 15. Plot selected GO terms for Syncytium-M2
###############################################################################

go_show_file <- file.path(go_dir, "Syncytium_M2_GO_show.txt")

if (file.exists(go_show_file)) {
  go_show <- read.delim(
    go_show_file,
    stringsAsFactors = FALSE
  )
} else {
  go_show <- Msyn_GO %>%
    filter(!is.na(p.adjust)) %>%
    group_by(Category) %>%
    arrange(p.adjust, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup()
}

if (nrow(go_show) > 0) {
  
  go_show <- go_show[order(go_show$p.adjust, decreasing = TRUE), ]
  go_show$term <- factor(go_show$Description, levels = go_show$Description)
  
  p_go <- ggplot(
    go_show,
    aes(term, -log10(p.adjust))
  ) +
    geom_col(aes(fill = Category), width = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = c("#3490DE", "#F07B3F", "#EA5455")) +
    facet_grid(Category ~ ., scales = "free_y", space = "free_y") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(color = "black", fill = "transparent")
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
    coord_flip() +
    labs(
      x = "",
      y = "-Log10 adjusted P value"
    )
  
  ggsave(
    filename = file.path(go_dir, "Syncytium_M2_GO_show.pdf"),
    plot = p_go,
    width = 10,
    height = 8
  )
}

###############################################################################
## 16. Helper functions for TOM-based gene network visualization
###############################################################################

format_gene_id <- function(gene_id) {
  ifelse(
    grepl("^ann1\\.", gene_id),
    gene_id,
    paste0("ann1.", gene_id)
  )
}

export_gene_tom_weights <- function(
    TOM_matrix,
    target_gene,
    output_file
) {
  target_gene <- format_gene_id(target_gene)
  
  if (!target_gene %in% rownames(TOM_matrix)) {
    warning("Target gene was not found in TOM: ", target_gene)
    return(NULL)
  }
  
  gene_weight <- as.data.frame(TOM_matrix[target_gene, , drop = FALSE])
  gene_weight <- as.data.frame(t(gene_weight))
  colnames(gene_weight) <- "weight"
  gene_weight$geneid <- rownames(gene_weight)
  
  gene_weight <- gene_weight %>%
    arrange(desc(weight))
  
  write.table(
    gene_weight,
    file = output_file,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  return(gene_weight)
}

get_network_genes <- function(
    TOM_matrix,
    target_genes,
    manual_gene_file = NULL,
    top_n = 30
) {
  target_genes <- format_gene_id(target_genes)
  
  if (!is.null(manual_gene_file) && file.exists(manual_gene_file)) {
    manual_genes <- read.table(
      manual_gene_file,
      header = FALSE,
      stringsAsFactors = FALSE
    )
    
    network_genes <- c(target_genes, manual_genes[, 1])
    network_genes <- format_gene_id(network_genes)
    network_genes <- unique(network_genes)
    network_genes <- intersect(network_genes, rownames(TOM_matrix))
    
    return(network_genes)
  }
  
  candidate_genes <- c()
  
  for (target_gene in target_genes) {
    if (!target_gene %in% rownames(TOM_matrix)) {
      next
    }
    
    weights <- TOM_matrix[target_gene, ]
    weights <- sort(weights, decreasing = TRUE)
    weights <- weights[names(weights) != target_gene]
    
    candidate_genes <- c(
      candidate_genes,
      names(weights)[seq_len(min(top_n, length(weights)))]
    )
  }
  
  network_genes <- unique(c(target_genes, candidate_genes))
  network_genes <- intersect(network_genes, rownames(TOM_matrix))
  
  return(network_genes)
}

plot_tom_network <- function(
    TOM_matrix,
    network_genes,
    output_file,
    width = 6.5,
    height = 5.5
) {
  network_genes <- intersect(network_genes, rownames(TOM_matrix))
  
  if (length(network_genes) < 2) {
    warning("Fewer than two genes were available for network plotting.")
    return(NULL)
  }
  
  cur_TOM <- TOM_matrix[network_genes, network_genes]
  
  graph <- cur_TOM %>%
    igraph::graph_from_adjacency_matrix(
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    ) %>%
    tidygraph::as_tbl_graph(directed = FALSE) %>%
    tidygraph::activate(nodes)
  
  p <- ggraph(graph, layout = layout_in_circle(graph)) +
    geom_edge_link(aes(alpha = weight), color = "#cccccc") +
    geom_node_point(color = "#0b6b9d", size = 3) +
    geom_node_label(
      aes(label = sub("^ann1\\.", "", name)),
      repel = TRUE,
      max.overlaps = Inf,
      fontface = "italic"
    ) +
    theme_void()
  
  ggsave(
    filename = output_file,
    plot = p,
    width = width,
    height = height
  )
  
  return(p)
}

###############################################################################
## 17. TOM-based network visualization for candidate genes
###############################################################################

candidate_genes <- c(
  "Glyma.04G100400",
  "Glyma.10G212900",
  "Glyma.15G013900"
)

candidate_gene_info <- data.frame(
  target_gene = candidate_genes,
  manual_file = file.path(
    manual_list_dir,
    c(
      "modTOM_Glyma.04G100400_net.txt",
      "modTOM_Glyma.10G212900_net.txt",
      "modTOM_Glyma.15G013900_net.txt"
    )
  ),
  output_prefix = c(
    "Glyma.04G100400",
    "Glyma.10G212900",
    "Glyma.15G013900"
  ),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(candidate_gene_info))) {
  
  target_gene <- candidate_gene_info$target_gene[i]
  output_prefix <- candidate_gene_info$output_prefix[i]
  manual_file <- candidate_gene_info$manual_file[i]
  
  target_gene_full <- format_gene_id(target_gene)
  
  weight_file <- file.path(
    network_dir,
    paste0("modTOM_", output_prefix, "_all.txt")
  )
  
  export_gene_tom_weights(
    TOM_matrix = TOM,
    target_gene = target_gene_full,
    output_file = weight_file
  )
  
  network_genes <- get_network_genes(
    TOM_matrix = TOM,
    target_genes = target_gene_full,
    manual_gene_file = manual_file,
    top_n = 30
  )
  
  write.table(
    data.frame(gene = sub("^ann1\\.", "", network_genes)),
    file = file.path(network_dir, paste0(output_prefix, "_network_genes_used.txt")),
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  )
  
  plot_tom_network(
    TOM_matrix = TOM,
    network_genes = network_genes,
    output_file = file.path(network_dir, paste0(output_prefix, "_related_net.pdf")),
    width = 6.5,
    height = 5.5
  )
}

###############################################################################
## 18. Joint network visualization for three candidate genes
###############################################################################

three_gene_manual_file <- file.path(
  manual_list_dir,
  "modTOM_three_genes_net.txt"
)

three_gene_network_genes <- get_network_genes(
  TOM_matrix = TOM,
  target_genes = candidate_genes,
  manual_gene_file = three_gene_manual_file,
  top_n = 20
)

write.table(
  data.frame(gene = sub("^ann1\\.", "", three_gene_network_genes)),
  file = file.path(network_dir, "three_genes_network_genes_used.txt"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

plot_tom_network(
  TOM_matrix = TOM,
  network_genes = three_gene_network_genes,
  output_file = file.path(network_dir, "three_genes_related_net.pdf"),
  width = 6.5,
  height = 5.5
)

###############################################################################
## 19. Save hdWGCNA object and workspace
###############################################################################

saveRDS(
  F3_syn,
  file = file.path(output_dir, "F3_syncytium_hdWGCNA_object.rds")
)

save.image(
  file = file.path(output_dir, "hdWGCNA_workflow.RData")
)