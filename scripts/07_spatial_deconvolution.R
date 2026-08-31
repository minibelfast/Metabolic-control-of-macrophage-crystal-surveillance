# ---------------------------------------------------------------------------
# 07_spatial_deconvolution.R
#
# Reference based deconvolution of the Visium spots with RCTD (spacexr), using
# the annotated scRNA-seq atlas of the same kidneys as reference.
#   * doublet mode -> dominant cell type per spot  (SpatialDimPlot)
#   * full mode    -> cell type weights per spot   (STdeconvolve scatterpies)
#
# The reference cell type panel differs between the two sections (see
# rctd_mode_by_group in R/config.R): the GLY section is deconvolved against the
# whole kidney panel, the NC section against the immune panel only.
#
# Input : data/processed/spatial_<group>_clustered.rda   (script 05)
#         data/scRNA/IntegrateData_3.Rdata
# Output: results/07_spatial_deconvolution/
#         data/processed/spatial_<group>_annotated.rda
#
# Run from the repository root:  Rscript scripts/07_spatial_deconvolution.R
# ---------------------------------------------------------------------------

source("R/setup.R")

suppressPackageStartupMessages({
  library(spacexr)
  library(STdeconvolve)
})

outdir <- file.path(dir_results, "07_spatial_deconvolution")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Build the single cell reference ------------------------------------

#' RCTD reference restricted to a set of cell types
#'
#' @param object    annotated scRNA-seq Seurat object (labels in cell_type)
#' @param cell_types cell types to keep, in the order used for plotting
build_reference <- function(object, cell_types) {
  Idents(object) <- "cell_type"
  object <- subset(object, idents = cell_types)
  object$cell_type <- factor(object$cell_type, levels = cell_types)

  counts   <- GetAssayData(object, slot = "counts", assay = "RNA")
  celltype <- droplevels(as.factor(object$cell_type))
  names(celltype) <- colnames(object)

  Reference(counts, celltype)
}

sc_all <- load_seurat(path_sc_all)

## ---- 2. Deconvolve each section --------------------------------------------

for (group in sample_groups) {

  mode_name  <- rctd_mode_by_group[[group]]
  cell_types <- rctd_reference_sets[[mode_name]]
  message("=== RCTD on ", group, " with the '", mode_name, "' reference ===")

  stRNA <- load_seurat(spatial_object_path(group, "clustered"), name = "stRNA")

  # Spot level counts and pixel coordinates define the RCTD query
  st_count <- GetAssayData(stRNA, slot = "counts", assay = "Spatial")
  st_loc   <- stRNA@images$spatial@coordinates[, c("row", "col")]
  query    <- SpatialRNA(st_loc, st_count)

  reference <- build_reference(sc_all, cell_types)
  RCTD      <- create.RCTD(query, reference, max_cores = 1)

  ## ---- 2a. doublet mode: dominant cell type per spot ---------------------
  # doublet : 1-2 cell types per spot, for high resolution assays
  # full    : any number of cell types per spot, recommended for Visium
  # multi   : extension of doublet mode allowing more than two cell types
  rctd_doublet <- run.RCTD(RCTD, doublet_mode = "doublet")
  results_df   <- rctd_doublet@results$results_df

  # Keep only the spots RCTD could assign, then transfer the label
  stRNA <- stRNA[, rownames(results_df)]
  stRNA$celltype1 <- results_df$first_type
  Idents(stRNA) <- "celltype1"

  save_plot(
    SpatialDimPlot(stRNA),
    file.path(outdir, sprintf("spatial_dimplot_%s_%s.pdf", group, mode_name)),
    width = 12, height = 12
  )
  write.csv(
    results_df,
    file.path(outdir, sprintf("rctd_doublet_%s_%s.csv", group, mode_name)),
    row.names = TRUE
  )

  # This object carries the deconvolved labels used by script 08
  save(stRNA, file = spatial_object_path(group, "annotated"))

  ## ---- 2b. full mode: cell type weights per spot -------------------------
  rctd_full   <- run.RCTD(RCTD, doublet_mode = "full")
  norm_weights <- normalize_weights(rctd_full@results$weights)
  write.csv(
    as.matrix(norm_weights),
    file.path(outdir, sprintf("rctd_full_weights_%s_%s.csv", group, mode_name)),
    row.names = TRUE
  )

  # Scatterpie representation of the weight matrix; STdeconvolve expects the
  # coordinates as columns named y / x.
  theta <- as.matrix(norm_weights)
  pos   <- st_loc[rownames(theta), ]
  colnames(pos) <- c("y", "x")

  p <- vizAllTopics(
    theta      = theta,
    pos        = pos,
    topicOrder = seq_len(ncol(theta)),
    topicCols  = rainbow(ncol(theta)),
    groups     = NA,
    group_cols = NA,
    r          = 0.6,   # scatterpie radius, depends on the pixel coordinates
    lwd        = 0.3,
    showLegend = TRUE,
    plotTitle  = sprintf("%s (%s reference)", group, mode_name)
  ) +
    guides(fill = guide_legend(ncol = 2))

  save_plot(
    p,
    file.path(outdir, sprintf("scatterpies_%s_%s.pdf", group, mode_name)),
    width = 12, height = 12
  )
}
