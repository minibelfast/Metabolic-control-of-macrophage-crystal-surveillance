# ---------------------------------------------------------------------------
# 06_spatial_feature_maps.R
#
# Spatial expression of the gene programmes discussed in the paper:
#   * fatty-acid uptake / transport
#   * crystal adhesion machinery (Cd81, integrins, Rac1, Arp2/3)
#   * ligand-receptor genes of the SPP1 / MIF / CCL axes
#   * metabolic markers
#   * a module score of the core crystal adhesion gene set
#
# Input : data/processed/spatial_<group>_clustered.rda   (script 05)
# Output: results/06_spatial_features/<group>/
#
# Run from the repository root:  Rscript scripts/06_spatial_feature_maps.R
# ---------------------------------------------------------------------------

source("R/setup.R")

outdir <- file.path(dir_results, "06_spatial_features")

for (group in sample_groups) {

  message("=== spatial feature maps: ", group, " ===")
  stRNA <- load_seurat(spatial_object_path(group, "clustered"), name = "stRNA")
  group_dir <- file.path(outdir, group)

  ## ---- 1. Fatty-acid transport programme ---------------------------------
  # A shared colour scale (keep.scale = "all") makes the sections comparable,
  # image.alpha = 0 hides the H&E background.
  save_spatial_features(
    stRNA, genes_fatty_acid_transport, group_dir,
    prefix = "fatty_acid", keep.scale = "all", image.alpha = 0
  )

  ## ---- 2. Metabolic markers ----------------------------------------------
  save_spatial_features(
    stRNA, genes_metabolic_markers, group_dir,
    prefix = "metabolic", keep.scale = "all", image.alpha = 0
  )

  ## ---- 3. Crystal adhesion machinery -------------------------------------
  save_spatial_features(
    stRNA, genes_crystal_adhesion, group_dir,
    prefix = "crystal_adhesion"
  )

  ## ---- 4. Ligand-receptor genes of the enriched axes ----------------------
  save_spatial_features(stRNA, genes_spp1_axis, group_dir, prefix = "spp1_axis")
  save_spatial_features(stRNA, genes_mif_axis,  group_dir, prefix = "mif_axis")
  save_spatial_features(stRNA, genes_ccl_axis,  group_dir, prefix = "ccl_axis")
  save_spatial_features(stRNA, genes_other_lr,  group_dir, prefix = "other_lr")

  ## ---- 5. Crystal adhesion module score ----------------------------------
  module_name <- names(geneset_crystal_module)[1]
  stRNA <- AddModuleScore(stRNA, features = geneset_crystal_module, name = module_name)
  # AddModuleScore appends the index of the gene set to the column name
  score_col <- paste0(module_name, "1")

  save_plot(
    SpatialFeaturePlot(stRNA, features = score_col),
    file.path(group_dir, sprintf("module_score_%s.pdf", module_name)),
    width = 12, height = 12
  )
}
