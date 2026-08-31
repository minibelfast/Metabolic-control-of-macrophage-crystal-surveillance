# ---------------------------------------------------------------------------
# 01_scrna_kidney_overview.R
#
# Whole-kidney scRNA-seq: marker based validation of the cell type labels,
# cell type composition between control (NC) and crystal nephropathy (GLY)
# kidneys, and expression of a single gene across all compartments.
#
# Input : path_sc_all  (annotated Seurat object, labels in meta.data$cell_type)
# Output: results/01_scrna_kidney/
# ---------------------------------------------------------------------------

source("R/setup.R")

outdir <- file.path(dir_results, "01_scrna_kidney")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sc_all <- load_seurat(path_sc_all)
DefaultAssay(sc_all) <- "RNA"

print(table(sc_all$cell_type, sc_all$orig.ident))

## ---- 1. Marker panel per compartment -------------------------------------

p_markers <- marker_violin_plot(
  sc_all,
  genes       = markers_kidney,
  cell_levels = kidney_cell_levels,
  group_by    = "cell_type"
)
save_plot(p_markers, file.path(outdir, "marker_violin.pdf"),
          width = 8, height = 15)

## ---- 2. Cell type composition --------------------------------------------

# Contribution of each condition to a given compartment
p_condition <- composition_bar_plot(
  sc_all,
  group_by      = "cell_type",
  sample_by     = "orig.ident",
  cell_levels   = kidney_cell_levels,
  sample_levels = group_levels,
  by            = "celltype"
)
save_plot(p_condition, file.path(outdir, "composition_by_celltype.pdf"),
          width = 6, height = 4)

# Composition of each condition
p_sample <- composition_bar_plot(
  sc_all,
  group_by      = "cell_type",
  sample_by     = "orig.ident",
  cell_levels   = kidney_cell_levels,
  sample_levels = group_levels,
  by            = "sample"
)
save_plot(p_sample, file.path(outdir, "composition_by_sample.pdf"),
          width = 6, height = 2)

## ---- 3. Expression of a single gene across compartments ------------------

# Macrophages versus proximal tubule cells, the two compartments that dominate
# the crystal response.
vln <- group_violin_plot(
  sc_all,
  gene         = "Triobp",
  group_by     = "cell_type",
  group_levels = kidney_cell_levels,
  comparisons  = list(c("Macro", "PT")),
  jitter_size  = 0.03
)
save_plot(vln$plot, file.path(outdir, "Triobp_by_celltype.pdf"),
          width = 6, height = 8)
write.csv(vln$data, file.path(outdir, "Triobp_by_celltype.csv"))

# Same gene inside the endothelial compartment only, NC versus GLY
Idents(sc_all) <- "cell_type"
endo <- subset(sc_all, idents = "Endo")

vln_endo <- group_violin_plot(
  endo,
  gene         = "Triobp",
  group_by     = "orig.ident",
  group_levels = group_levels,
  comparisons  = list(group_levels),
  jitter_size  = 3,
  drop_zero    = TRUE
)
save_plot(vln_endo$plot, file.path(outdir, "Triobp_endothelium_NC_vs_GLY.pdf"),
          width = 4, height = 4)
write.csv(vln_endo$data,
          file.path(outdir, "Triobp_endothelium_NC_vs_GLY.csv"))
