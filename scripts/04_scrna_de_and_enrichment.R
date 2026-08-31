# ---------------------------------------------------------------------------
# 04_scrna_de_and_enrichment.R
#
# Differential expression in the myeloid compartment and functional annotation
# of the genes that mark Cx3cr1+ macrophages:
#   (a) Cx3cr1+ macrophages versus the other macrophage subsets
#   (b) Cx3cr1+ macrophages, crystal nephropathy (GLY) versus control (NC)
#   (c) GO biological process enrichment of the genes up in Cx3cr1+ macrophages
#
# Input : path_sc_macro
# Output: results/04_scrna_de_enrichment/
# ---------------------------------------------------------------------------

source("R/setup.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(stringr)
})

outdir <- file.path(dir_results, "04_scrna_de_enrichment")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sc_macro <- load_seurat(path_sc_macro)
DefaultAssay(sc_macro) <- "RNA"

#' One-versus-one differential expression with a signed fold change
#'
#' FindAllMarkers() is run with only.pos = TRUE, so every gene is reported for
#' the group it is enriched in. The sign of log2FoldChange is flipped for the
#' reference group, which turns the two one-sided tables into a single table
#' that can be shown on a volcano plot.
#'
#' @param object    Seurat object, already restricted to the cells of interest
#' @param group_col meta.data column holding the two groups
#' @param positive  group that keeps the positive fold change
#' @return data.frame with the columns expected by plot_volcano()
signed_markers <- function(object, group_col, positive) {
  Idents(object) <- group_col
  markers <- FindAllMarkers(
    object,
    test.use        = "wilcox",
    only.pos        = TRUE,
    logfc.threshold = de_logfc_threshold,
    min.pct         = de_min_pct
  )
  markers$log2FoldChange <- ifelse(markers$cluster == positive,
                                   markers$avg_log2FC, -markers$avg_log2FC)
  markers$padj <- markers$p_val_adj
  markers
}

## ---- 1. Cx3cr1+ versus the other macrophage subsets ----------------------

Idents(sc_macro) <- "celltype"
macro <- subset(sc_macro,
                idents = c("Cx3cr1+ Macro", "Spp1+ Macro", "Birc5+ Macro"))
macro$Group <- ifelse(macro$celltype == "Cx3cr1+ Macro", "Cx3cr1+", "Cx3cr1-")

markers_subset <- signed_markers(macro, "Group", positive = "Cx3cr1+")
write.csv(markers_subset,
          file.path(outdir, "markers_Cx3cr1_vs_other_macro.csv"),
          row.names = FALSE)

save_plot(plot_volcano(markers_subset,
                       logfc = volcano_logfc, padj = volcano_padj,
                       label_genes = volcano_label_genes),
          file.path(outdir, "volcano_Cx3cr1_vs_other_macro.pdf"),
          width = 5, height = 4)

## ---- 2. Cx3cr1+ macrophages, GLY versus NC -------------------------------

cx3cr1 <- subset(sc_macro, idents = "Cx3cr1+ Macro")

markers_condition <- signed_markers(cx3cr1, "orig.ident", positive = "NC")
write.csv(markers_condition,
          file.path(outdir, "markers_Cx3cr1_NC_vs_GLY.csv"),
          row.names = FALSE)

save_plot(plot_volcano(markers_condition,
                       logfc = volcano_logfc, padj = volcano_padj,
                       label_genes = volcano_label_genes),
          file.path(outdir, "volcano_Cx3cr1_NC_vs_GLY.pdf"),
          width = 5, height = 4)

## ---- 3. GO enrichment of the Cx3cr1+ signature ---------------------------

signature <- markers_subset$gene[
  markers_subset$log2FoldChange > volcano_logfc &
    markers_subset$padj < volcano_padj
]

entrez <- unlist(mget(signature, org.Mm.egSYMBOL2EG, ifnotfound = NA))
entrez <- unique(entrez[!is.na(entrez)])

go <- enrichGO(entrez, OrgDb = org.Mm.eg.db, ont = "BP",
               pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
go_df <- as.data.frame(go)
go_df <- go_df[go_df$pvalue < 0.05 & go_df$qvalue < 0.05, ]
write.table(go_df, file.path(outdir, "GO_BP.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Terms highlighted in the manuscript figure
go_selected <- go_df[go_df$ID %in% go_terms_of_interest, ]
save_plot(go_bubble_plot(go_selected),
          file.path(outdir, "GO_BP_selected.pdf"), width = 9, height = 4)

# Table formatted for GOplot (chord diagrams etc.)
go_plot_input <- go_selected[, c("ID", "Description", "p.adjust", "geneID")]
names(go_plot_input) <- c("ID", "Term", "adj_pval", "Genes")
go_plot_input$Genes    <- str_replace_all(go_plot_input$Genes, "/", ",")
go_plot_input$Category <- "BP"
write.csv(go_plot_input, file.path(outdir, "GO_BP_selected_goplot.csv"),
          row.names = FALSE)
