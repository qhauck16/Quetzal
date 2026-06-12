#!/usr/bin/env Rscript
# Quetzal v1.0 -- per-gene DE analysis.
#
# (res.RDS, matrix.RDS) -> de_res.RDS via fastTopics::de_analysis on the
# poisson2multinom-converted fit. Same call as the DE block in v0.1
# gene_plots_and_objs.R; factored out here so make_plots.R is the only
# place that paints, and so this step doesn't run at all in genome_wide
# mode (DE results aren't used for the genome-wide flashier collapse).

# null-coalesce helper (Rscript on R < 4.4 doesn't ship one).
`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(optparse)
})
if (!requireNamespace("fastTopics", quietly = TRUE)) {
  remotes::install_github("stephenslab/fastTopics",
                           upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(fastTopics)
})

option_list <- list(
  make_option("--res",    type = "character",
              help = "res.RDS (fastTopics Poisson NMF fit from fit_pnmf.R)"),
  make_option("--matrix", type = "character",
              help = "matrix.RDS (canonical counts)"),
  make_option("--output", type = "character",
              help = "de_res.RDS output"),
  make_option("--threads", type = "integer", default = 1L,
              help = "control nc passed to de_analysis [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$res), !is.null(opt$matrix), !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

mat <- readRDS(opt$matrix)
res <- readRDS(opt$res)

# Pass-through if either upstream was skipped.
if (isTRUE(mat$skipped) || isTRUE(res$skipped)) {
  reason <- if (isTRUE(mat$skipped)) mat$reason else res$reason
  gene   <- mat$gene_name %||% res$gene_name
  saveRDS(list(skipped = TRUE, reason = reason, gene_name = gene),
           opt$output)
  message(sprintf("  [pass-through %s] upstream skipped (%s)", gene, reason))
  quit(save = "no", status = 0)
}

# fastTopics::de_analysis wants the poisson2multinom-converted fit and
# the same counts matrix used to fit it (samples x junctions).
multi <- poisson2multinom(res)
dat   <- t(as.matrix(mat$counts))    # samples x junctions

de_res <- de_analysis(multi, dat, control = list(nc = opt$threads))
de_res$gene_name <- mat$gene_name

saveRDS(de_res, opt$output)
message(sprintf("  [ok %s] DE done; wrote %s", mat$gene_name, opt$output))
