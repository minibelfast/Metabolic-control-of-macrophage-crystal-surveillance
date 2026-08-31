# ---------------------------------------------------------------------------
# config.R
# Central configuration: input/output paths, sample groups and gene sets.
# All analysis scripts are meant to be executed from the repository root, so
# every path below is relative to that root.
# ---------------------------------------------------------------------------

## ---- Directories ---------------------------------------------------------

# Raw / processed input data (not tracked by git, see README for the layout)
dir_data <- "data"

# scRNA-seq Seurat objects
# - sc_all   : all kidney cell types, cell labels stored in meta.data$cell_type
# - sc_macro : re-clustered myeloid/macrophage compartment, labels in meta.data$celltype
path_sc_all   <- file.path(dir_data, "scRNA", "IntegrateData_3.Rdata")
path_sc_macro <- file.path(dir_data, "scRNA", "IntegrateData_3_Macro.Rdata")

# Spatial transcriptomics (10x Visium) input directories, one per group.
# Each directory must contain filtered_feature_bc_matrix.h5 and spatial/
dir_spatial <- file.path(dir_data, "spatial")

# External Visium validation sample used for the metabolic analysis (GSE206306)
dir_spatial_validation <- file.path(dir_data, "spatial_validation", "GSM6250307")

# Intermediate Seurat objects written by the pipeline
dir_processed <- file.path(dir_data, "processed")
dir.create(dir_processed, showWarnings = FALSE, recursive = TRUE)

#' Path of an intermediate spatial object
#'
#' @param group "NC" or "GLY"
#' @param stage "clustered" (script 05) or "annotated" (script 07)
spatial_object_path <- function(group, stage = c("clustered", "annotated")) {
  stage <- match.arg(stage)
  file.path(dir_processed, sprintf("spatial_%s_%s.rda", group, stage))
}

# Output directory for figures and tables
dir_results <- "results"
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

## ---- Experimental design -------------------------------------------------

# NC  : control mice
# GLY : glyoxylate-induced calcium oxalate crystal nephropathy
sample_groups <- c("NC", "GLY")

# Order used consistently on all plots (control first)
group_levels <- c("NC", "GLY")

## ---- Cell type ordering --------------------------------------------------

# Whole-kidney compartments (meta.data$cell_type of the sc_all object)
kidney_cell_levels <- c(
  "T cell", "B cell", "Neutro", "Macro", "CD-IC", "CD-PC",
  "DCT", "LOH", "PT", "Podo", "Endo", "unknown"
)

# Myeloid subsets (meta.data$celltype of the sc_macro object)
myeloid_cell_levels <- c(
  "Itgae+ DC", "Ccl17+ DC", "Cd8a+ T cell", "Klrb1c+ T cell",
  "Birc5+ Macro", "Spp1+ Macro", "Spp1+ Mono", "Cd163+ DC",
  "Stard9+ Mono", "Cx3cr1+ Macro"
)

# Cell types used as RCTD reference. Two deconvolution modes are reported:
# "tissue" resolves the whole kidney, "immune" focuses on infiltrating cells.
rctd_reference_sets <- list(
  tissue = c("T cell", "CD-IC", "CD-PC", "DCT", "LOH",
             "Macro", "Neutro", "Podo", "PT", "B cell"),
  immune = c("T cell", "Macro", "Neutro", "B cell")
)

# Reference used for each section in the manuscript figures
rctd_mode_by_group <- c(GLY = "tissue", NC = "immune")

## ---- Gene sets -----------------------------------------------------------

# Marker panel for the whole-kidney annotation (violin plot)
markers_kidney <- c(
  "Kdr", "Nrp1", "Nphs1", "Nphs2", "Slc27a2", "Lrp2", "Slc12a1", "Umod",
  "Slc12a3", "Pvalb", "Aqp2", "Hsd11b2", "Atp6v1g3", "Atp6v0d2", "Insrr",
  "Rhbg", "Itgam", "Adgre1", "S100a8", "S100a9", "Cd79a", "Cd79b", "Ltb",
  "Cxcr6", "Nkg7"
)

# Marker panel for the myeloid subsets (bubble plot)
markers_myeloid <- c(
  "Irf8", "Itgae", "Ly75", "Pdcd1lg2", "Ccl17", "Gzma", "Cd8a", "Cd4",
  "Cd3g", "Cd3e", "Klrb1c", "Cd68", "Hmgb2", "Birc5", "Ly6c1", "Il4ra",
  "Itgam", "Spp1", "Umod", "Slc12a3", "Cd163", "Macf1", "Gm26917",
  "Stard9", "C1qa", "C1qb", "Cd81"
)

# Fatty-acid uptake / transport programme
genes_fatty_acid_transport <- c(
  "Cd36",
  paste0("Fabp", 1:7),
  paste0("Slc27a", 1:6)
)

# Crystal adhesion / phagocytic cup machinery
genes_crystal_adhesion <- c("Cd81", "Itgb1", "Itgb5", "Rac1", "Actr2", "Actr3")

# Module score gene set (core crystal adhesion module)
geneset_crystal_module <- list(crystal_adhesion = c("Cd81", "Itgb1", "Itgb5", "Rac1"))

# Ligand / receptor genes of the spatially enriched pathways
genes_spp1_axis <- c("Spp1", "Itgav", "Itgb6", "Itgb1", "Itga4")
genes_mif_axis  <- c("Mif", "Cd74", "Cd44", "Cxcr4", "Ackr3")
genes_ccl_axis  <- c("Ccl3", "Ccl4", "Ccl6", "Ccr2", "Ccr5")
genes_other_lr  <- c("App", "Col4a4", "Sdc4")

# Metabolic markers shown on the tissue sections
genes_metabolic_markers <- c("Pkm", "Nat10", "Miox", "Ifng")

# Genes highlighted on the volcano plots
volcano_label_genes <- c("Cd81", "C1qa", "C1qc", "C1qb", "Tmem176b")

# GO biological processes selected for the enrichment bubble plot
go_terms_of_interest <- c(
  "GO:1905517", "GO:0002495", "GO:0022409", "GO:0050870", "GO:0002720",
  "GO:0045785", "GO:0022604", "GO:0032760", "GO:1900271", "GO:0051017",
  "GO:2001280"
)

# CellChat signalling pathway highlighted in the spatial network plots and the
# ligand / receptor pair of that pathway shown as spot level expression
# CellChatDB.mouse can be restricted to one annotation category
# ("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"); NA uses all.
cellchat_db_subset           <- "Cell-Cell Contact"
cellchat_pathway_of_interest <- "JAM"
cellchat_genes_of_interest   <- c("Timp1", "Jam3")

# Ligand-receptor pairs visualised in space (CellChat naming)
lr_pairs_of_interest <- c(
  "SPP1_ITGAV_ITGB6", "SPP1_ITGAV_ITGB5", "SPP1_ITGAV_ITGB1",
  "SPP1_ITGA4_ITGB1",
  "MIF_CD74_CXCR4", "MIF_CD74_CD44", "MIF_ACKR3",
  "CCL6_CCR2", "COL4A4_SDC4", "APP_CD74"
)

# KEGG pathway scored with scMetabolism on the external validation section
metabolic_pathway_of_interest <- "Drug metabolism - cytochrome P450"

## ---- Analysis parameters -------------------------------------------------

# Number of principal components used downstream (chosen with choose_pcs())
n_pcs_scrna   <- 24
n_pcs_spatial <- 19

# Differential expression thresholds
de_logfc_threshold <- 0.25
de_min_pct         <- 0.1

# Thresholds for volcano plots and enrichment input
volcano_logfc <- 1
volcano_padj  <- 0.05
