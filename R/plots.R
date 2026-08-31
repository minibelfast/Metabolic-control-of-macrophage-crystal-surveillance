# ---------------------------------------------------------------------------
# plots.R
# Re-usable figure builders. Each function returns a ggplot object so that the
# analysis scripts only decide *what* to plot and *where* to save it.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(reshape2)
  library(RColorBrewer)
  library(ggsignif)
})

#' Marker gene bubble plot
#'
#' Scaled expression is averaged per cell type and combined with the fraction of
#' cells expressing the gene (dot size).
#'
#' @param object    Seurat object
#' @param genes     character vector of genes (y axis, plotted top to bottom)
#' @param cell_levels ordering of the cell types on the x axis
#' @param group_by  meta.data column holding the cell labels
#' @param assay     assay to take normalised expression from
#' @param cap       upper bound applied to the mean scaled expression
#' @return ggplot object
marker_bubble_plot <- function(object, genes, cell_levels,
                               group_by = "celltype", assay = "RNA", cap = 3) {
  expr <- Seurat::GetAssayData(object, slot = "data", assay = assay)[genes, ]
  scaled <- as.data.frame(scale(t(as.matrix(expr))))
  scaled[[group_by]] <- object@meta.data[rownames(scaled), group_by]

  plot_df <- scaled %>%
    tidyr::pivot_longer(cols = all_of(genes), names_to = "gene", values_to = "expr") %>%
    group_by(.data[[group_by]], gene) %>%
    summarise(
      exp   = mean(expr),
      ratio = sum(expr > min(expr)) / length(expr),
      .groups = "drop"
    ) %>%
    rename(celltype = all_of(group_by))

  plot_df$celltype <- factor(plot_df$celltype, levels = rev(cell_levels))
  plot_df$gene     <- factor(plot_df$gene, levels = rev(genes))
  plot_df$exp      <- pmin(plot_df$exp, cap)

  ggplot(plot_df, aes(celltype, gene, size = ratio, colour = exp)) +
    geom_point() +
    scale_x_discrete("") +
    scale_y_discrete("") +
    scale_colour_gradientn(
      colours = rev(c("#FFD92F", "#FEE391", brewer.pal(11, "Spectral")[7:11]))
    ) +
    scale_size_continuous(limits = c(0, 1)) +
    theme_classic() +
    theme(axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 90))
}

#' Stacked marker gene violin plot (one facet per gene)
#'
#' @param object Seurat object
#' @param genes  character vector of genes (facets, top to bottom)
#' @param cell_levels ordering of the cell types on the x axis
#' @param group_by meta.data column holding the cell labels
#' @param assay  assay to take normalised expression from
#' @return ggplot object
marker_violin_plot <- function(object, genes, cell_levels,
                               group_by = "cell_type", assay = "RNA") {
  expr <- as.data.frame(Seurat::GetAssayData(object, slot = "data", assay = assay)[genes, ])
  expr$gene <- rownames(expr)

  plot_df <- melt(expr, id.vars = "gene",
                  variable.name = "barcode", value.name = "expr")
  plot_df$group <- object@meta.data[as.character(plot_df$barcode), group_by]
  plot_df$group <- factor(plot_df$group, levels = cell_levels)
  plot_df$gene  <- factor(plot_df$gene, levels = genes)

  ggplot(plot_df, aes(group, expr)) +
    geom_violin(aes(fill = group), scale = "width") +
    facet_grid(gene ~ ., scales = "free_y") +
    scale_fill_brewer(palette = "Set3") +
    scale_x_discrete("") +
    scale_y_continuous("") +
    theme_bw() +
    theme(
      axis.text.x.bottom = element_text(angle = 0, hjust = 0.5, vjust = 0.5),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none"
    )
}

#' Stacked bar plot of cell type composition
#'
#' @param object Seurat object
#' @param group_by meta.data column with the cell labels
#' @param sample_by meta.data column with the sample / condition
#' @param cell_levels ordering of the cell labels
#' @param sample_levels ordering of the samples
#' @param by one of "celltype" (samples split within each cell type) or
#'   "sample" (cell types split within each sample)
#' @param palette colour vector
#' @return ggplot object
composition_bar_plot <- function(object, group_by, sample_by,
                                 cell_levels, sample_levels,
                                 by = c("celltype", "sample"),
                                 palette = qualitative_palette()) {
  by <- match.arg(by)

  plot_df <- melt(table(object@meta.data[[group_by]],
                        object@meta.data[[sample_by]]))
  colnames(plot_df) <- c("Cluster", "Sample", "Number")
  plot_df$Cluster <- factor(plot_df$Cluster, levels = cell_levels)
  plot_df$Sample  <- factor(plot_df$Sample, levels = sample_levels)

  mapping <- if (by == "celltype") {
    aes(x = Number, y = Cluster, fill = Sample)
  } else {
    aes(x = Number, y = Sample, fill = Cluster)
  }
  n_fill <- if (by == "celltype") length(sample_levels) else length(cell_levels)

  ggplot(plot_df, mapping) +
    geom_bar(stat = "identity", width = 0.8, position = "fill") +
    scale_fill_manual(values = palette[seq_len(n_fill)]) +
    labs(x = "", y = "Ratio") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(size = 12, colour = "black"),
      axis.text.x = element_text(size = 12, colour = "black"),
      axis.text.x.bottom = element_text(hjust = 1, vjust = 1, angle = 45)
    )
}

#' Single gene UMAP feature plot with a grey-to-purple gradient
#'
#' @param object Seurat object with a "umap" reduction
#' @param gene   gene symbol
#' @param assay  assay to take normalised expression from
#' @return ggplot object
umap_feature_plot <- function(object, gene, assay = "RNA") {
  plot_df <- as.data.frame(Embeddings(object, "umap"))
  plot_df$expr <- fetch_expression(object, gene, assay)[rownames(plot_df)]

  ggplot(plot_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = expr)) +
    scale_colour_gradient(low = "grey", high = "purple") +
    scale_x_continuous("") +
    scale_y_continuous("") +
    ggtitle(gene) +
    theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      legend.position = "left",
      plot.title = element_text(hjust = 0.5, size = 14)
    )
}

#' Violin plot of one gene across groups, with pairwise significance testing
#'
#' Replaces the three near-identical violin plot blocks of the original
#' notebooks: the grouping variable, the compared levels and the jitter size are
#' now arguments.
#'
#' @param object Seurat object (subset beforehand if needed)
#' @param gene   gene symbol
#' @param group_by meta.data column used on the x axis
#' @param group_levels ordering of the groups
#' @param comparisons list of length-2 character vectors passed to geom_signif()
#' @param jitter_size size of the jittered points
#' @param drop_zero drop cells with zero expression before plotting
#' @param assay  assay to take normalised expression from
#' @return list with the ggplot object (`plot`) and the plotted data (`data`)
group_violin_plot <- function(object, gene, group_by, group_levels,
                              comparisons, jitter_size = 0.03,
                              drop_zero = FALSE, assay = "RNA") {
  expr <- fetch_expression(object, gene, assay)
  plot_df <- object@meta.data
  plot_df$expr <- expr[rownames(plot_df)]
  if (drop_zero) plot_df <- plot_df[plot_df$expr > 0, ]
  plot_df$group <- factor(plot_df[[group_by]], levels = group_levels)
  plot_df <- plot_df[!is.na(plot_df$group), ]

  p <- ggplot(plot_df, aes(group, expr, fill = group)) +
    geom_violin() +
    geom_signif(comparisons = comparisons, step_increase = 0.1,
                map_signif_level = FALSE, test = t.test,
                size = 1, textsize = 6) +
    geom_jitter(shape = 16, size = jitter_size,
                position = position_jitter(0.2)) +
    labs(x = "Group", y = "Expression", title = gene) +
    guides(fill = guide_legend(title = "Group")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
      axis.text.x = element_text(size = 15, face = "bold",
                                 angle = 25, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = 18, face = "bold"),
      axis.title.x = element_text(size = 0, face = "bold"),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.line.x = element_line(linetype = 1, colour = "black", size = 1),
      axis.line.y = element_line(linetype = 1, colour = "black", size = 1),
      legend.text = element_text(size = 15)
    )

  list(plot = p, data = plot_df)
}

#' GO enrichment bubble plot
#'
#' @param go_df data.frame from clusterProfiler with Description, Count, pvalue
#' @return ggplot object
go_bubble_plot <- function(go_df) {
  ggplot(go_df, aes(Count, Description)) +
    geom_point(aes(size = Count, colour = -log10(pvalue))) +
    cols4all::scale_color_continuous_c4a_seq("rd_pu") +
    labs(x = "Number of genes", y = "") +
    theme_bw() +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 11),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 11),
      plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
    )
}

#' Write one spatial feature plot per gene
#'
#' Replaces the dozens of copy-pasted SpatialFeaturePlot / ggsave blocks of the
#' original notebook.
#'
#' @param object  spatial Seurat object
#' @param genes   character vector of genes; genes absent from the object are skipped
#' @param outdir  directory the pdf files are written to
#' @param prefix  file name prefix, e.g. "NC_fatty_acid"
#' @param ...     further arguments forwarded to SpatialFeaturePlot()
#' @return character vector of the files written
save_spatial_features <- function(object, genes, outdir, prefix, ...) {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  present <- intersect(genes, rownames(object))
  missing <- setdiff(genes, present)
  if (length(missing) > 0) {
    message("skipping genes absent from the object: ",
            paste(missing, collapse = ", "))
  }

  vapply(present, function(gene) {
    p <- SpatialFeaturePlot(object, features = gene, ...)
    file <- file.path(outdir, sprintf("%s_%s.pdf", prefix, gene))
    ggsave(file, plot = p, width = 12, height = 12, bg = "white")
    file
  }, character(1))
}
