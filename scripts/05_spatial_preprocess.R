# ---------------------------------------------------------------------------
# 05_spatial_preprocess.R
#
# Visium spatial transcriptomics: normalisation, dimensionality reduction and
# unsupervised clustering of the control (NC) and crystal nephropathy (GLY)
# kidney sections.
#
# Input : data/spatial/<group>/filtered_feature_bc_matrix.h5 + spatial/
# Output: results/05_spatial_preprocess/  (PC diagnostics, cluster maps)
#         data/processed/spatial_<group>_clustered.rda
#
# Run from the repository root:  Rscript scripts/05_spatial_preprocess.R
# ---------------------------------------------------------------------------

source("R/setup.R")

outdir <- file.path(dir_results, "05_spatial_preprocess")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

for (group in sample_groups) {

  message("=== processing spatial sample: ", group, " ===")
  data_dir <- file.path(dir_spatial, group)

  ## ---- 1. Load the Visium section ----------------------------------------
  stRNA <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    slice    = "spatial"
  )

  ## ---- 2. Normalisation (SCTransform) and PCA ----------------------------
  stRNA <- SCTransform(stRNA, assay = "Spatial", verbose = FALSE)
  stRNA <- RunPCA(stRNA, assay = "SCT", verbose = FALSE)

  # Diagnostic only: the number of PCs actually used is fixed in config.R so
  # that the analysis stays reproducible.
  pcs <- choose_pcs(stRNA, plot_file = file.path(outdir, sprintf("pc_selection_%s.pdf", group)))
  message("suggested number of PCs: ", pcs,
          " (config.R uses n_pcs_spatial = ", n_pcs_spatial, ")")

  ## ---- 3. Clustering and UMAP --------------------------------------------
  dims <- seq_len(n_pcs_spatial)
  stRNA <- FindNeighbors(stRNA, reduction = "pca", dims = dims)
  stRNA <- FindClusters(stRNA, verbose = FALSE)
  stRNA <- RunUMAP(stRNA, reduction = "pca", dims = dims)

  ## ---- 4. Cluster maps ---------------------------------------------------
  save_plot(
    SpatialPlot(stRNA, label = TRUE, label.size = 5),
    file.path(outdir, sprintf("spatial_clusters_%s.pdf", group)),
    width = 12, height = 12
  )
  save_plot(
    DimPlot(stRNA, reduction = "umap", label = TRUE),
    file.path(outdir, sprintf("umap_clusters_%s.pdf", group)),
    width = 8, height = 7
  )

  ## ---- 5. Store the clustered object for the downstream scripts ----------
  save(stRNA, file = spatial_object_path(group, "clustered"))
}
