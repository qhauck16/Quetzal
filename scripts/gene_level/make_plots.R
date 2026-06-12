#!/usr/bin/env Rscript
# Quetzal v1.0 -- per-gene whole_factor.html (and optional de_factor.html).
#
# Ports the plotting block of v0.1 gene_plots_and_objs.R (lines 303-544)
# behind a configurable-grouping CLI. Inputs:
#   --res, --de_res    : RDS outputs of fit_pnmf.R and de_analysis.R
#   --matrix           : matrix.RDS (for junction_info -> per-junction coords)
#   --gencode          : .RDS exon model (gene rendered via ggbio::geom_alignment)
#   --gene_name        : gene to plot
#   --whole_output     : path to whole_factor.html
#   --de_output        : path to de_factor.html (only written when at least one
#                         factor has DE-significant junctions; otherwise the
#                         file is a tiny placeholder so Snakemake's output
#                         declaration is satisfied)
#
# Optional grouping for the structure_plot:
#   --sample_metadata     : TSV with one row per sample
#   --sample_id_column    : column in metadata holding sample IDs that
#                            match rownames(res$L)            [rail_id]
#   --grouping_column     : column in metadata used to group samples
#   When --grouping_column is unset OR --sample_metadata is unset, the
#   structure plot draws one ungrouped block (matches v0.1 when run on
#   non-TCGA data).
#
# DE thresholds (used to gate the de_factor.html junction arcs):
#   --lfc_thresh, --abs_thresh, --sval_thresh

# null-coalesce helper.
`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(ggbio)
  library(patchwork)
  library(plotly)
  library(htmlwidgets)
  library(purrr)
  library(GenomicRanges)
})
if (!requireNamespace("fastTopics", quietly = TRUE)) {
  remotes::install_github("stephenslab/fastTopics",
                           upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(fastTopics)
})

option_list <- list(
  make_option("--res",          type = "character"),
  make_option("--de_res",       type = "character"),
  make_option("--matrix",       type = "character"),
  make_option("--gencode",      type = "character",
              default = "data/hg38_granges.RDS"),
  make_option("--gene_name",    type = "character"),
  make_option("--whole_output", type = "character"),
  make_option("--de_output",    type = "character"),
  make_option("--sample_metadata",  type = "character", default = NA_character_),
  make_option("--sample_id_column", type = "character", default = "rail_id"),
  make_option("--grouping_column",  type = "character", default = NA_character_),
  make_option("--lfc_thresh",   type = "double", default = 2),
  make_option("--abs_thresh",   type = "double", default = 0.01),
  make_option("--sval_thresh",  type = "double", default = 0.05)
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$res), !is.null(opt$matrix),
           !is.null(opt$whole_output), !is.null(opt$de_output),
           !is.null(opt$gene_name))
dir.create(dirname(opt$whole_output), showWarnings = FALSE, recursive = TRUE)
dir.create(dirname(opt$de_output),    showWarnings = FALSE, recursive = TRUE)

placeholder <- function(path, reason) {
  writeLines(sprintf("<!-- %s skipped: %s -->", opt$gene_name, reason), path)
}

# ---- load inputs ---------------------------------------------------------

mat    <- readRDS(opt$matrix)
res    <- readRDS(opt$res)
de_res <- readRDS(opt$de_res)

if (isTRUE(mat$skipped) || isTRUE(res$skipped) || isTRUE(de_res$skipped)) {
  reason <- (if (isTRUE(mat$skipped)) mat$reason
             else if (isTRUE(res$skipped)) res$reason
             else de_res$reason)
  placeholder(opt$whole_output, sprintf("upstream skipped (%s)", reason))
  placeholder(opt$de_output,    sprintf("upstream skipped (%s)", reason))
  message(sprintf("  [pass-through %s] %s", opt$gene_name, reason))
  quit(save = "no", status = 0)
}

multi       <- poisson2multinom(res)
counts      <- mat$counts
junction_info <- mat$junction_info
samples     <- colnames(counts)
num_factors <- ncol(multi$L)

# ---- gene model from gencode --------------------------------------------

gencode   <- readRDS(opt$gencode)
gene_data <- gencode[gencode$gene_name == opt$gene_name &
                        gencode$type     == "exon"]
if (length(gene_data) == 0L) {
  placeholder(opt$whole_output, "gene not in gencode")
  placeholder(opt$de_output,    "gene not in gencode")
  message(sprintf("  [skip %s] gene not in gencode", opt$gene_name))
  quit(save = "no", status = 0)
}
t_models <- rtracklayer::split(gene_data, gene_data$transcript_id)

round_down_to_1000 <- function(x) floor(x / 1000) * 1000
round_up_to_1000   <- function(x) ceiling(x / 1000) * 1000
min_coord <- round_down_to_1000(min(start(gene_data), end(gene_data)))
max_coord <- round_up_to_1000(max(start(gene_data), end(gene_data)))

base_p <- ggplot() +
  ggbio::geom_alignment(t_models, gap.geom = "arrow", label = FALSE,
                         exon.rect.h = 0.4) +
  theme_minimal() +
  scale_x_continuous(breaks = seq(min_coord, max_coord, by = 10000))

build <- ggplot_build(base_p)
max_y <- max(build$layout$panel_params[[1]]$y.range)

# ---- optional grouping for structure_plot --------------------------------

grouping_vec <- NULL
if (!is.na(opt$grouping_column) && !is.na(opt$sample_metadata) &&
     nzchar(opt$grouping_column) && nzchar(opt$sample_metadata)) {
  meta <- read_tsv(opt$sample_metadata, show_col_types = FALSE)
  if (!opt$sample_id_column %in% colnames(meta)) {
    message(sprintf("  WARN: sample metadata missing column '%s'; skipping grouping",
                     opt$sample_id_column))
  } else if (!opt$grouping_column %in% colnames(meta)) {
    message(sprintf("  WARN: sample metadata missing column '%s'; skipping grouping",
                     opt$grouping_column))
  } else {
    sid_to_group <- setNames(as.character(meta[[opt$grouping_column]]),
                              as.character(meta[[opt$sample_id_column]]))
    grouping_vec <- sid_to_group[samples]
    if (anyNA(grouping_vec)) {
      grouping_vec[is.na(grouping_vec)] <- "(unmapped)"
    }
  }
}

# ---- palette ------------------------------------------------------------

palet <- c("red", "blue", "green", "purple", "orange", "cyan2", "goldenrod4",
            "pink", "magenta", "grey", "brown", "darkolivegreen", "skyblue2",
            "yellow", "chartreuse4", "black",
            # extras if num_factors goes past 16
            "tomato", "navy", "limegreen", "violet")
palet <- palet[seq_len(max(num_factors, 1L))]
palet_labels <- setNames(palet, as.character(seq_len(length(palet))))

# ---- structure plot -----------------------------------------------------

struc_plot <- structure_plot(multi,
                              grouping = grouping_vec,
                              gap      = 100,
                              colors   = palet[seq_len(num_factors)],
                              n        = 5000) +
  theme(aspect.ratio = 1 / 5) +
  labs(title = paste0("Factorisation of GENE: ", opt$gene_name))

# ---- factor x junction layer dataframes ---------------------------------

# Per-junction (factor, loading, start, end). multi$F is junctions x factors;
# rownames(F) is the canonical "chr:start-end:strand" set by fit_pnmf.R.
# But the coords are also in mat$junction_info, which is cleaner.
F_df <- as.data.frame(multi$F)
colnames(F_df) <- as.character(seq_len(num_factors))
F_df$junction_id <- rownames(multi$F)
F_df <- F_df %>%
  inner_join(junction_info %>% dplyr::select(junction_id, start, end),
              by = "junction_id")

s_vals_mat <- de_res$svalue
s_vals_mat[is.na(s_vals_mat)] <- 1
lfc_mat    <- de_res$postmean
lfc_mat[is.na(lfc_mat)] <- 0

whole_line_df <- data.frame(x = numeric(), y = numeric(),
                             xend = numeric(), yend = numeric(),
                             width = numeric(), factor = character())
line_df       <- whole_line_df

k <- 1
for (j in seq_len(num_factors)) {
  width_col <- F_df[[as.character(j)]]
  this_whole <- data.frame(x = F_df$start, y = max_y + 5 * j,
                            xend = F_df$end, yend = max_y + 5 * j,
                            width = width_col,
                            factor = as.character(j)) %>%
    filter(width > opt$abs_thresh / 10)
  whole_line_df <- rbind(whole_line_df, this_whole)

  max_lfc <- apply(lfc_mat, 1, max)
  s_for_j <- s_vals_mat[, j]
  l_for_j <- lfc_mat[,    j]
  # de-significant: width > abs, sval < threshold, lfc > thresh and is the
  # row's max lfc (i.e. the junction prefers this factor)
  keep_de <- (width_col > opt$abs_thresh) &
              (s_for_j < opt$sval_thresh) &
              (l_for_j > opt$lfc_thresh) &
              (l_for_j == max_lfc)
  if (any(keep_de)) {
    this_de <- data.frame(x = F_df$start[keep_de], y = max_y + 5 * k,
                           xend = F_df$end[keep_de], yend = max_y + 5 * k,
                           width = width_col[keep_de],
                           factor = as.character(j))
    line_df <- rbind(line_df, this_de)
    k <- k + 1
  }
}

# ---- whole factor plot --------------------------------------------------

running_whole_plot <- base_p
height_scaling     <- 0.7
if (nrow(whole_line_df) > 0L) {
  for (i in seq_len(nrow(whole_line_df))) {
    y_base   <- max_y + as.numeric(whole_line_df$factor[i]) * 4
    x_vals   <- seq(whole_line_df$x[i], whole_line_df$xend[i],
                     length.out = 100)
    span     <- whole_line_df$xend[i] - whole_line_df$x[i]
    arc      <- sin((x_vals - whole_line_df$x[i]) / span * pi) * height_scaling
    y_top    <- arc + y_base +
                  whole_line_df$width[i] / max(whole_line_df$width) * 3.5
    y_bottom <- arc + y_base
    polygon_df <- data.frame(x = c(x_vals, rev(x_vals)),
                              y = c(y_top,  rev(y_bottom)))
    running_whole_plot <- running_whole_plot +
      geom_polygon(data = polygon_df, aes(x = x, y = y),
                    fill = palet[as.numeric(whole_line_df$factor[i])],
                    alpha = 1)
  }
}
running_whole_plot <- running_whole_plot +
  scale_y_continuous(breaks = c(0, 10), labels = c("", ""))

ggp_whole_factor <- ggplotly(running_whole_plot, dynamicTicks = TRUE)
if (nrow(whole_line_df) > 0L) {
  whole_line_df$label <- paste0(whole_line_df$x, "-", whole_line_df$xend)
  polygon_traces <- 4:(3 + nrow(whole_line_df))
  ggp_whole_factor <- purrr::reduce(
    seq_along(polygon_traces),
    function(plotly_obj, i) {
      plotly::style(plotly_obj, hoverinfo = "text",
                     text = whole_line_df$label[i],
                     traces = polygon_traces[i])
    },
    .init = ggp_whole_factor
  )
}

# Strand arrow
strand_char <- as.character(t_models@unlistData@strand[1])
if (strand_char == "+") {
  arrow_start <- min_coord
  arrow_end   <- min_coord + (max_coord - min_coord) / 10
} else {
  arrow_start <- max_coord
  arrow_end   <- max_coord - (max_coord - min_coord) / 10
}

ggp_struc <- ggplotly(struc_plot)

ggp_whole_factor <- ggp_whole_factor %>%
  add_annotations(x = arrow_end, y = -1.4, ax = arrow_start, ay = -1.4,
                   xref = "x2", yref = "y2", axref = "x2", ayref = "y2",
                   text = "", showarrow = TRUE, arrowhead = 10,
                   arrowsize = 3, arrowwidth = 0.5, arrowcolor = "black") %>%
  add_annotations(x = (arrow_start + arrow_end) / 2, y = -0.3,
                   text = "gene direction", showarrow = FALSE,
                   font = list(size = 9), xanchor = "center")

htmlwidgets::saveWidget(
  subplot(ggp_struc, ggp_whole_factor, nrows = 2, margin = 0.07),
  opt$whole_output)
message(sprintf("  [ok %s] wrote %s", opt$gene_name, opt$whole_output))

# ---- DE factor plot (optional) ------------------------------------------

if (nrow(line_df) == 0L) {
  placeholder(opt$de_output, "no DE-significant junctions")
  message(sprintf("  [skip de %s] no DE-significant junctions", opt$gene_name))
  unlink(sprintf("%s_files", tools::file_path_sans_ext(opt$whole_output)),
          recursive = TRUE)
  quit(save = "no", status = 0)
}

running_de_plot <- base_p
for (i in seq_len(nrow(line_df))) {
  y_base   <- max_y + as.numeric(line_df$factor[i]) * 4
  x_vals   <- seq(line_df$x[i], line_df$xend[i], length.out = 100)
  span     <- line_df$xend[i] - line_df$x[i]
  arc      <- sin((x_vals - line_df$x[i]) / span * pi) * height_scaling
  y_top    <- arc + y_base +
                line_df$width[i] / max(whole_line_df$width) * 3.5
  y_bottom <- arc + y_base
  polygon_df <- data.frame(x = c(x_vals, rev(x_vals)),
                            y = c(y_top,  rev(y_bottom)))
  running_de_plot <- running_de_plot +
    geom_polygon(data = polygon_df, aes(x = x, y = y),
                  fill = palet[as.numeric(line_df$factor[i])], alpha = 1)
}
running_de_plot <- running_de_plot +
  scale_y_continuous(breaks = c(0, 10), labels = c("", ""))

ggp_de_factor <- ggplotly(running_de_plot, dynamicTicks = TRUE)
line_df$label <- paste0(line_df$x, "-", line_df$xend)
polygon_traces <- 4:(3 + nrow(line_df))
ggp_de_factor <- purrr::reduce(
  seq_along(polygon_traces),
  function(plotly_obj, i) {
    plotly::style(plotly_obj, hoverinfo = "text",
                   text = line_df$label[i],
                   traces = polygon_traces[i])
  },
  .init = ggp_de_factor
)
ggp_de_factor <- ggp_de_factor %>%
  add_annotations(x = arrow_end, y = -1.4, ax = arrow_start, ay = -1.4,
                   xref = "x2", yref = "y2", axref = "x2", ayref = "y2",
                   text = "", showarrow = TRUE, arrowhead = 10,
                   arrowsize = 3, arrowwidth = 0.5, arrowcolor = "black") %>%
  add_annotations(x = (arrow_start + arrow_end) / 2, y = -0.3,
                   text = "gene direction", showarrow = FALSE,
                   font = list(size = 9), xanchor = "center")

htmlwidgets::saveWidget(
  subplot(ggp_struc, ggp_de_factor, nrows = 2, margin = 0.07),
  opt$de_output)
message(sprintf("  [ok %s] wrote %s", opt$gene_name, opt$de_output))

# clean up htmlwidgets sidecar dirs
unlink(sprintf("%s_files", tools::file_path_sans_ext(opt$whole_output)),
        recursive = TRUE)
unlink(sprintf("%s_files", tools::file_path_sans_ext(opt$de_output)),
        recursive = TRUE)
