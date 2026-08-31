# ---------------------------------------------------------------------------
# utils.R
# Small helpers shared by the single-cell and spatial pipelines.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(RColorBrewer)
})

#' Build a large qualitative colour palette
#'
#' Concatenates every qualitative RColorBrewer palette, which yields ~73
#' reasonably distinguishable colours - enough for all clusters used here.
#'
#' @return character vector of hex colours
qualitative_palette <- function() {
  info <- brewer.pal.info[brewer.pal.info$category == "qual", ]
  unlist(mapply(brewer.pal, info$maxcolors, rownames(info)))
}

#' Load a Seurat object stored in an .Rdata file
#'
#' The archives shipped with this project store the object under the name
#' "scRNA1". Loading into a temporary environment keeps the global workspace
#' clean and makes the dependency of every script explicit.
#'
#' @param path path to the .Rdata file
#' @param name name of the object inside the archive
#' @return Seurat object
load_seurat <- function(path, name = "scRNA1") {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!name %in% ls(env)) {
    stop(sprintf("'%s' does not contain an object called '%s' (found: %s)",
                 path, name, paste(ls(env), collapse = ", ")))
  }
  get(name, envir = env)
}

#' Choose the number of principal components to keep
#'
#' Applies three criteria and keeps the most conservative one:
#'   1. cumulative variance explained > 90%
#'   2. variance explained by the PC itself < 5%
#'   3. difference between two consecutive PCs < 0.1%
#'
#' @param object   Seurat object with a "pca" reduction
#' @param plot_file optional path; if given the elbow diagnostic is written there
#' @return integer, the suggested number of PCs
choose_pcs <- function(object, plot_file = NULL) {
  pct  <- object[["pca"]]@stdev / sum(object[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)

  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:(length(pct) - 1)] - pct[2:length(pct)]) > 0.1),
              decreasing = TRUE)[1] + 1
  pcs <- min(co1, co2)

  if (!is.null(plot_file)) {
    plot_df <- data.frame(pct = pct, cumu = cumu, rank = seq_along(pct))
    p <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) +
      geom_text() +
      geom_vline(xintercept = 90, color = "grey") +
      geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      theme_bw()
    ggsave(plot_file, plot = p, width = 8, height = 8)
  }

  pcs
}

#' Save a ggplot object with the defaults used throughout the project
#'
#' @param plot   ggplot object
#' @param file   output path (pdf)
#' @param width,height figure size in inches
save_plot <- function(plot, file, width = 8, height = 7) {
  ggsave(file, plot = plot, width = width, height = height, bg = "white")
  invisible(file)
}

#' Fetch the expression vector of a single gene from the RNA assay
#'
#' Works with both Seurat v4 slots and v5 layers.
#'
#' @param object Seurat object
#' @param gene   gene symbol
#' @param assay  assay name
#' @return named numeric vector of normalised expression
fetch_expression <- function(object, gene, assay = "RNA") {
  mat <- Seurat::GetAssayData(object, slot = "data", assay = assay)
  if (!gene %in% rownames(mat)) {
    stop(sprintf("gene '%s' is absent from assay '%s'", gene, assay))
  }
  mat[gene, ]
}
