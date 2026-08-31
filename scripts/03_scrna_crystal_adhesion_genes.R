# ---------------------------------------------------------------------------
# 03_scrna_crystal_adhesion_genes.R
#
# Expression of the crystal adhesion machinery in the myeloid compartment:
# UMAP feature plots of the individual genes and, within Cx3cr1+ macrophages,
# the induction of Itgb5 in crystal nephropathy.
#
# Input : path_sc_macro
# Output: results/03_scrna_crystal_genes/
# ---------------------------------------------------------------------------

source("R/setup.R")

outdir <- file.path(dir_results, "03_scrna_crystal_genes")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sc_macro <- load_seurat(path_sc_macro)
DefaultAssay(sc_macro) <- "RNA"

if (!"umap" %in% Reductions(sc_macro)) {
  sc_macro <- RunUMAP(sc_macro, dims = seq_len(n_pcs_scrna))
}

## ---- 1. UMAP feature plots ------------------------------------------------

for (gene in genes_crystal_adhesion) {
  if (!gene %in% rownames(sc_macro)) {
    message("skipping gene absent from the object: ", gene)
    next
  }
  save_plot(umap_feature_plot(sc_macro, gene),
            file.path(outdir, sprintf("umap_%s.pdf", gene)),
            width = 8, height = 7)
}

## ---- 2. Itgb5 in Cx3cr1+ macrophages, NC versus GLY ----------------------

Idents(sc_macro) <- "celltype"
cx3cr1 <- subset(sc_macro, idents = "Cx3cr1+ Macro")

# Only expressing cells are kept, so that the comparison reports the level of
# the gene rather than the fraction of positive cells.
vln <- group_violin_plot(
  cx3cr1,
  gene         = "Itgb5",
  group_by     = "orig.ident",
  group_levels = group_levels,
  comparisons  = list(group_levels),
  jitter_size  = 3,
  drop_zero    = TRUE
)
save_plot(vln$plot, file.path(outdir, "Itgb5_Cx3cr1_macro_NC_vs_GLY.pdf"),
          width = 4, height = 4)
write.csv(vln$data, file.path(outdir, "Itgb5_Cx3cr1_macro_NC_vs_GLY.csv"))
