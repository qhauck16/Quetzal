#!/usr/bin/env Rscript
# Quetzal v1.0 -- per-gene unannotated-factor extraction.
#
# Mirrors lines 116-138 of v0.1 `fasttopics_to_flashier.R` but as a
# Snakemake rule per gene so the work parallelises across genes / nodes
# instead of running in one serial R loop.
#
# Inputs (both produced upstream by canonicalise + fit_pnmf):
#   --matrix : matrix.RDS  -- supplies junction_info$annotated (0/1) per junction
#   --res    : res.RDS     -- the fastTopics Poisson NMF fit
#
# Output:
#   --output : unann_factors.tsv -- wide TSV with one row per sample and
#              one column per kept factor, named "<gene>.k<N>" using the
#              factor's column index in res$L. Always includes a `sample_id`
#              column. When the gene was skipped upstream OR no factors
#              pass the loading-ratio threshold, the TSV holds only the
#              `sample_id` column (no factor columns), so the genome_wide
#              `softimpute_flash` gather step can full_join indiscriminately
#              and the gene contributes nothing.
#
# Factor-selection rule (carried verbatim from v0.1):
#   keep factor k  iff  max(res$F[unann, k]) * loading_ratio
#                          > max(res$F[, k])
# i.e. the factor's top-loading unannotated junction must exceed
# (1 / loading_ratio) of the factor's top overall loading.
# Default loading_ratio = 2 means "top unannotated >= 50% of top overall".

# null-coalesce helper (Rscript on R < 4.4 doesn't ship one).
`%||%` <- function(a, b) if (!is.null(a)) a else b

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
})
# fastTopics::poisson2multinom is what makes res$L sum to 1 per sample
# (needed because the genome_wide flashier collapse expects a multinomial L).
if (!requireNamespace("fastTopics", quietly = TRUE)) {
  remotes::install_github("stephenslab/fastTopics",
                           upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(fastTopics)
})

option_list <- list(
  make_option("--matrix", type = "character",
              help = "matrix.RDS (canonical per-gene; carries annotated 0/1)"),
  make_option("--res",    type = "character",
              help = "res.RDS (fastTopics fit from fit_pnmf.R)"),
  make_option("--output", type = "character",
              help = "unann_factors.tsv output path"),
  make_option("--unann_factor_loading_ratio", type = "double", default = 2.0,
              help = "keep factor k iff max(F[unann, k]) * ratio > max(F[, k]) [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$matrix), !is.null(opt$res), !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

mat <- readRDS(opt$matrix)
res <- readRDS(opt$res)

# Helper: write a sample_id-only TSV (factor columns absent so a downstream
# full_join contributes no new columns for this gene).
write_empty <- function(samples, reason) {
  tibble(sample_id = samples) %>% write_tsv(opt$output)
  message(sprintf("  [empty %s] %s",
                   if (!is.null(mat$gene_name)) mat$gene_name else "?",
                   reason))
}

# --- pass-through skips ----------------------------------------------------

if (isTRUE(mat$skipped)) {
  write_empty(character(0), sprintf("upstream matrix skipped: %s", mat$reason))
  quit(save = "no", status = 0)
}
if (isTRUE(res$skipped)) {
  write_empty(rownames(mat$counts) %||% character(0),
              sprintf("upstream res skipped: %s", res$reason))
  quit(save = "no", status = 0)
}

# --- locate unannotated junctions + filter F ------------------------------

annotated_status <- mat$junction_info$annotated
if (anyNA(annotated_status)) {
  # classify_junctions wasn't run; refuse rather than guess.
  write_empty(rownames(res$L) %||% character(0),
              "junction_info$annotated has NA values; run classify_junctions.R first")
  quit(save = "no", status = 0)
}

unann_ids <- mat$junction_info$junction_id[annotated_status == 0L]
if (length(unann_ids) == 0L) {
  write_empty(rownames(res$L), "no unannotated junctions in this gene")
  quit(save = "no", status = 0)
}

# fit_pnmf set rownames(res$F) = junction_id, so this is a clean lookup.
mask     <- rownames(res$F) %in% unann_ids
if (!any(mask)) {
  write_empty(rownames(res$L),
              "no junction_id overlap between matrix and res (unexpected)")
  quit(save = "no", status = 0)
}
unann_F  <- res$F[mask, , drop = FALSE]

# --- factor selection (v0.1 rule) -----------------------------------------

max_unann   <- apply(unann_F,  2, max)
max_overall <- apply(res$F,    2, max)
keep        <- (max_unann * opt$unann_factor_loading_ratio) > max_overall

if (!any(keep)) {
  write_empty(rownames(res$L), "no factor passed unann_factor_loading_ratio")
  quit(save = "no", status = 0)
}

# --- multinomial L -> kept-factor loadings --------------------------------

multi <- poisson2multinom(res)
L     <- as.data.frame(multi$L[, keep, drop = FALSE])
gene  <- res$gene_name %||% mat$gene_name
colnames(L) <- paste0(gene, ".k", which(keep))

out <- L %>%
  rownames_to_column("sample_id") %>%
  as_tibble()

write_tsv(out, opt$output)
message(sprintf("  [ok %s] kept %d / %d factors (%s); wrote %s",
                 gene, sum(keep), length(keep),
                 paste(which(keep), collapse = ","),
                 opt$output))
