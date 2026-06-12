#!/usr/bin/env Rscript
# Quetzal v1.0 -- gather per-gene unann_factors.tsv files and fit the
# genome-wide softImpute + flashier matrix factorisation.
#
# Replaces the back half of v0.1 fasttopics_to_flashier.R (lines 200-217),
# minus the TCGA-specific QC block (RIN / %C / avgQ / unique-mapped % /
# productive-unproductive). v1.0 keeps no QC -- pre-filter your inputs
# or post-filter the flashier object yourself.
#
# Inputs:
#   --input_dir : root of the genome_wide output dir
#                  (the Snakefile passes config['output_dir']/genome_wide).
#                  Walked recursively for any file named exactly
#                  'unann_factors.tsv', so empty / skipped per-gene files
#                  participate too (they just contribute no new feature
#                  columns at the full_join step).
#   --output    : softimpute_flash.RDS path.
#
# Knobs:
#   --sample_fraction : minimum fraction of TOTAL samples (union across all
#                       per-gene TSVs) a gene must cover to contribute its
#                       factors to the matrix. Default 0.8.
#   --greedy_kmax     : flashier::flash greedy_Kmax. Default 300.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(softImpute)
})
# flashier isn't always reliably packaged for conda-forge; install once
# into the active conda env on first use.
if (!requireNamespace("flashier", quietly = TRUE)) {
  remotes::install_github("willwerscheid/flashier",
                           upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(flashier)
})

option_list <- list(
  make_option("--input_dir",       type = "character",
              help = "root dir containing chr*/<gene>/unann_factors.tsv"),
  make_option("--output",          type = "character",
              help = "softimpute_flash.RDS output path"),
  make_option("--sample_fraction", type = "double",  default = 0.8,
              help = "min fraction of TOTAL samples a gene must cover [%default]"),
  make_option("--greedy_kmax",     type = "integer", default = 300L,
              help = "flashier greedy_Kmax [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input_dir), !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

# ---- gather per-gene TSVs ------------------------------------------------

tsv_paths <- list.files(opt$input_dir,
                         pattern    = "^unann_factors\\.tsv$",
                         recursive  = TRUE,
                         full.names = TRUE)
message(sprintf("Found %d per-gene unann_factors.tsv files under %s",
                 length(tsv_paths), opt$input_dir))
if (length(tsv_paths) == 0L) {
  stop("no unann_factors.tsv files found under --input_dir; did the upstream rules fire?")
}

per_gene <- lapply(tsv_paths, function(p) {
  df <- suppressWarnings(read_tsv(p, show_col_types = FALSE))
  # Empty stub: only sample_id column, no factor columns -> drop.
  if (!"sample_id" %in% colnames(df)) return(NULL)
  if (ncol(df) <= 1L)                 return(NULL)
  df
})
per_gene <- Filter(Negate(is.null), per_gene)
message(sprintf("  %d files carry at least one factor column",
                 length(per_gene)))

# ---- sample-coverage filter ----------------------------------------------

all_samples   <- unique(unlist(lapply(per_gene, function(d) d$sample_id)))
total_samples <- length(all_samples)
threshold     <- opt$sample_fraction * total_samples

n_before <- length(per_gene)
per_gene <- Filter(function(d) nrow(d) >= threshold, per_gene)
message(sprintf("  sample_fraction filter: %d / %d genes pass (>= %.2f * %d = %.0f samples covered)",
                 length(per_gene), n_before, opt$sample_fraction,
                 total_samples, threshold))
if (length(per_gene) == 0L) {
  stop("no gene passed the sample_fraction filter")
}

# ---- full_join into one sample x feature matrix --------------------------

big <- tibble(sample_id = all_samples)
for (d in per_gene) big <- full_join(big, d, by = "sample_id")

feature_df <- as.data.frame(big[, -1L, drop = FALSE])
rownames(feature_df) <- big$sample_id
X <- as.matrix(feature_df)
message(sprintf("  joined matrix: %d samples x %d features",
                 nrow(X), ncol(X)))

# ---- softImpute + flashier ----------------------------------------------

message("softImpute (rank.max = min(dim) - 1, lambda = 0) ...")
fit_si    <- softImpute(X, rank.max = min(dim(X)) - 1L, lambda = 0)
X_imputed <- complete(X, fit_si)
X_scaled  <- scale(X_imputed)

message(sprintf("flashier::flash (greedy_Kmax = %d, ebnm = point_exponential + point_laplace, var_type = 2, backfit = FALSE) ...",
                 opt$greedy_kmax))
fl <- flash(X_scaled,
             greedy_Kmax = opt$greedy_kmax,
             ebnm_fn     = list(ebnm_point_exponential, ebnm_point_laplace),
             var_type    = 2,
             backfit     = FALSE)
message(sprintf("flashier fit: %d factors", fl$n_factors))

saveRDS(fl, opt$output)
message(sprintf("wrote %s", opt$output))
