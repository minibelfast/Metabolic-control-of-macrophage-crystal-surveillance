# ---------------------------------------------------------------------------
# 02_scrna_myeloid_subsets.R
#
# Myeloid compartment of the kidney scRNA-seq data: choose the number of
# principal components, recompute the UMAP embedding and validate the ten
# myeloid subsets with a marker bubble plot.
#
# Input : path_sc_macro (myeloid Seurat object, labels in meta.data$celltype)
# Output: results/02_scrna_myeloid/
# ---------------------------------------------------------------------------

source("R/setup.R")

outdir <- file.path(dir_results, "02_scrna_myeloid")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sc_macro <- load_seurat(path_sc_macro)
DefaultAssay(sc_macro) <- "RNA"

print(table(sc_macro$celltype, sc_macro$orig.ident))

## ---- 1. Number of principal components ------------------------------------

# The diagnostic plot is written next to the figures; n_pcs_scrna in config.R
# records the value that was used for the manuscript.
pcs <- choose_pcs(sc_macro, plot_file = file.path(outdir, "pc_selection.pdf"))
message("suggested number of PCs: ", pcs,
        " (config.R uses n_pcs_scrna = ", n_pcs_scrna, ")")

## ---- 2. UMAP embedding ----------------------------------------------------

sc_macro <- RunUMAP(sc_macro, dims = seq_len(n_pcs_scrna))
write.csv(Embeddings(sc_macro, "umap"), file.path(outdir, "umap_embedding.csv"))

Idents(sc_macro) <- "celltype"
save_plot(DimPlot(sc_macro, reduction = "umap", label = TRUE),
          file.path(outdir, "umap_celltype.pdf"), width = 8, height = 7)

## ---- 3. Marker panel per subset ------------------------------------------

p_bubble <- marker_bubble_plot(
  sc_macro,
  genes       = markers_myeloid,
  cell_levels = myeloid_cell_levels,
  group_by    = "celltype"
)
save_plot(p_bubble, file.path(outdir, "marker_bubble.pdf"),
          width = 4.5, height = 12)

## ---- 4. Cell type composition of the myeloid compartment -----------------

p_sample <- composition_bar_plot(
  sc_macro,
  group_by      = "celltype",
  sample_by     = "orig.ident",
  cell_levels   = myeloid_cell_levels,
  sample_levels = group_levels,
  by            = "sample"
)
save_plot(p_sample, file.path(outdir, "composition_by_sample.pdf"),
          width = 6, height = 2)
