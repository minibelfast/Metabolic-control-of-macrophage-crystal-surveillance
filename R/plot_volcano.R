# ---------------------------------------------------------------------------
# plot_volcano.R
# Volcano plot used for the differential expression results.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
})

#' Volcano plot of a differential expression table
#'
#' @param result data.frame with the columns `gene`, `log2FoldChange` and `padj`
#' @param logfc  absolute log2 fold change cut-off (dashed vertical lines)
#' @param padj   adjusted p-value cut-off (dashed horizontal line)
#' @param label_genes optional character vector of genes to annotate
#' @return ggplot object
plot_volcano <- function(result, logfc = 1, padj = 0.05, label_genes = NULL) {
  stopifnot(all(c("gene", "log2FoldChange", "padj") %in% colnames(result)))

  result$change <- "NONE"
  result$change[result$log2FoldChange >  logfc & result$padj < padj] <- "UP"
  result$change[result$log2FoldChange < -logfc & result$padj < padj] <- "DOWN"

  colours <- c(NONE = "grey", UP = "red", DOWN = "blue")
  xlim <- max(abs(result$log2FoldChange), na.rm = TRUE)
  ylim <- max(-log10(result$padj), na.rm = TRUE)

  p <- ggplot(result, aes(log2FoldChange, -log10(padj))) +
    geom_point(aes(colour = change)) +
    geom_vline(xintercept = c(-logfc, logfc), lty = 2) +
    geom_hline(yintercept = -log10(padj), lty = 2) +
    scale_x_continuous(limits = c(-xlim, xlim)) +
    coord_fixed(ratio = (2 * xlim) / ylim) +
    scale_colour_manual(values = colours) +
    labs(x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(colour = "black"))

  if (!is.null(label_genes)) {
    p <- p +
      geom_label_repel(
        data = result[result$gene %in% label_genes, ],
        aes(label = gene, fill = change),
        colour = "white", fontface = "italic"
      ) +
      scale_fill_manual(values = colours)
  }

  p
}
