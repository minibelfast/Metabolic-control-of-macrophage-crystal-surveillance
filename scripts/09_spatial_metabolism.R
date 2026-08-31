# ---------------------------------------------------------------------------
# 09_spatial_metabolism.R
#
# Metabolic pathway activity of an independent Visium kidney section
# (GSE206306, sample GSM6250307) scored with scMetabolism (VISION / KEGG).
# This external section is used to confirm that the metabolic programme seen in
# our own sections is not specific to the GSE269465 dataset.
#
# Input : data/spatial_validation/GSM6250307/
# Output: results/09_spatial_metabolism/
#
# Run from the repository root:  Rscript scripts/09_spatial_metabolism.R
# ---------------------------------------------------------------------------

source("R/setup.R")

suppressPackageStartupMessages({
  library(scMetabolism)
})

outdir <- file.path(dir_results, "09_spatial_metabolism")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Load and cluster the validation section ----------------------------
stRNA <- Load10X_Spatial(
  data.dir = dir_spatial_validation,
  filename = "filtered_feature_bc_matrix.h5",
  slice    = "spatial"
)

stRNA <- SCTransform(stRNA, assay = "Spatial", verbose = FALSE)
stRNA <- RunPCA(stRNA, assay = "SCT", verbose = FALSE)
stRNA <- FindNeighbors(stRNA, reduction = "pca", dims = seq_len(n_pcs_spatial))
stRNA <- FindClusters(stRNA, verbose = FALSE)

save_plot(
  SpatialPlot(stRNA, label = TRUE, label.size = 5),
  file.path(outdir, "spatial_clusters.pdf"),
  width = 12, height = 12
)

## ---- 2. Metabolic pathway scoring ------------------------------------------
# scMetabolism reads the counts from an assay called "RNA", so the SCT assay is
# exposed under that name.
stRNA@assays$RNA <- stRNA@assays$SCT

stRNA <- sc.metabolism.Seurat(
  obj             = stRNA,
  method          = "VISION",
  imputation      = FALSE,
  ncores          = 2,
  metabolism.type = "KEGG"
)

scores <- stRNA@assays[["METABOLISM"]][["score"]]
writeLines(rownames(scores), file.path(outdir, "kegg_pathways.txt"))
write.csv(scores, file.path(outdir, "metabolism_scores.csv"), row.names = TRUE)

## ---- 3. Pathway of interest ------------------------------------------------
pathway <- metabolic_pathway_of_interest
stopifnot(pathway %in% rownames(scores))

save_plot(
  DotPlot.metabolism(obj = stRNA, pathway = pathway,
                     phenotype = "seurat_clusters", norm = "y"),
  file.path(outdir, "dotplot_by_cluster.pdf"),
  width = 6, height = 4
)

# Transfer the pathway score to the metadata so that it can be mapped in space
stRNA$pathway <- as.numeric(scores[pathway, ])
save_plot(
  SpatialFeaturePlot(stRNA, features = "pathway"),
  file.path(outdir, "spatial_pathway_score.pdf"),
  width = 12, height = 12
)
