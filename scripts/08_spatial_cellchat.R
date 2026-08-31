# ---------------------------------------------------------------------------
# 08_spatial_cellchat.R
#
# Spatially informed cell-cell communication analysis (CellChat) on the RCTD
# annotated Visium sections. The "Cell-Cell Contact" branch of CellChatDB.mouse
# is used, because the question is which contacts drive crystal surveillance.
#
# Input : data/processed/spatial_<group>_annotated.rda   (script 07)
#         data/spatial/<group>/spatial/scalefactors_json.json
# Output: results/08_spatial_cellchat/<group>/
#
# Run from the repository root:  Rscript scripts/08_spatial_cellchat.R
# ---------------------------------------------------------------------------

source("R/setup.R")

suppressPackageStartupMessages({
  library(CellChat)
  library(jsonlite)
})

outdir <- file.path(dir_results, "08_spatial_cellchat")

for (group in sample_groups) {

  message("=== CellChat on ", group, " ===")
  stRNA <- load_seurat(spatial_object_path(group, "annotated"), name = "stRNA")
  group_dir <- file.path(outdir, group)
  dir.create(group_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- 1. Build the spatial CellChat object -------------------------------
  expr <- GetAssayData(stRNA, slot = "data", assay = "SCT")
  meta <- data.frame(celltype = Idents(stRNA), row.names = names(Idents(stRNA)))

  spatial_locs <- GetTissueCoordinates(stRNA, scale = NULL,
                                       cols = c("imagerow", "imagecol"))

  visium_scale <- fromJSON(txt = file.path(dir_spatial, group, "spatial",
                                           "scalefactors_json.json"))
  # spot.diameter is the physical Visium spot size in micrometres; spot.diameter
  # and spot are required, the remaining entries are optional.
  scale_factors <- list(
    spot.diameter = 65,
    spot          = visium_scale$spot_diameter_fullres,
    fiducial      = visium_scale$fiducial_diameter_fullres,
    hires         = visium_scale$tissue_hires_scalef,
    lowres        = visium_scale$tissue_lowres_scalef
  )

  cellchat <- createCellChat(
    object        = expr,
    meta          = meta,
    group.by      = "celltype",
    datatype      = "spatial",
    coordinates   = spatial_locs,
    scale.factors = scale_factors
  )

  ## ---- 2. Ligand-receptor database ---------------------------------------
  db <- CellChatDB.mouse
  if (!is.na(cellchat_db_subset)) {
    db <- subsetDB(db, search = cellchat_db_subset, key = "annotation")
  }
  cellchat@DB <- db

  ## ---- 3. Preprocessing --------------------------------------------------
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  # Mouse data must be projected on the mouse PPI network
  cellchat <- projectData(cellchat, PPI.mouse)

  ## ---- 4. Communication probabilities ------------------------------------
  cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
  # Drop interactions supported by very few spots
  cellchat <- filterCommunication(cellchat, min.cells = 3)

  write.csv(subsetCommunication(cellchat),
            file.path(group_dir, "communication_lr.csv"), row.names = FALSE)
  write.csv(subsetCommunication(cellchat, slot.name = "netP"),
            file.path(group_dir, "communication_pathway.csv"), row.names = FALSE)

  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  ## ---- 5. Aggregated interaction network ---------------------------------
  # Node size is the number of spots per cell type, edge width the number of
  # ligand-receptor pairs (left) or the interaction strength (right).
  group_size <- as.numeric(table(cellchat@idents))
  pdf(file.path(group_dir, "network_circle.pdf"), width = 12, height = 6)
  par(mfrow = c(1, 2), xpd = TRUE)
  netVisual_circle(cellchat@net$count, vertex.weight = group_size,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Number of interactions")
  netVisual_circle(cellchat@net$weight, vertex.weight = group_size,
                   weight.scale = TRUE, label.edge = FALSE,
                   title.name = "Interaction weights/strength")
  dev.off()

  ## ---- 6. Pathway of interest --------------------------------------------
  pathways <- cellchat@netP$pathways
  writeLines(pathways, file.path(group_dir, "enriched_pathways.txt"))

  pathway <- cellchat_pathway_of_interest
  if (pathway %in% pathways) {

    pdf(file.path(group_dir, sprintf("network_spatial_%s.pdf", pathway)),
        width = 12, height = 12)
    netVisual_aggregate(cellchat, signaling = pathway, layout = "spatial",
                        edge.width.max = 2, vertex.size.max = 1,
                        alpha.image = 0.2, vertex.label.cex = 3.5)
    dev.off()

    save_plot(
      netAnalysis_contribution(cellchat, signaling = pathway),
      file.path(group_dir, sprintf("lr_contribution_%s.pdf", pathway)),
      width = 6, height = 4
    )

    present <- intersect(cellchat_genes_of_interest, rownames(expr))
    if (length(present) > 0) {
      save_plot(
        spatialFeaturePlot(cellchat, features = present, point.size = 0.8,
                           color.heatmap = "Reds", direction = 1),
        file.path(group_dir, sprintf("expression_%s.pdf", pathway)),
        width = 12, height = 6
      )
    }
  } else {
    message("pathway '", pathway, "' is not enriched in ", group, ", skipping")
  }

  ## ---- 7. Spatial maps of the selected ligand-receptor pairs -------------
  available_pairs <- rownames(cellchat@LR$LRsig)
  for (pair in lr_pairs_of_interest) {
    if (!pair %in% available_pairs) {
      message("ligand-receptor pair '", pair, "' is unavailable in ", group,
              ", skipping")
      next
    }
    save_plot(
      spatialFeaturePlot(cellchat, pairLR.use = pair, point.size = 3,
                         do.binary = TRUE, cutoff = 0.05,
                         enriched.only = FALSE, direction = 1),
      file.path(group_dir, sprintf("lr_pair_%s.pdf", pair)),
      width = 12, height = 12
    )
  }

  saveRDS(cellchat, file.path(dir_processed, sprintf("cellchat_%s.rds", group)))
}
