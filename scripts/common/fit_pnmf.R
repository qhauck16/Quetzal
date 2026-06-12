#!/usr/bin/env Rscript
# Quetzal v1.0 -- canonical matrix.RDS -> Poisson NMF res.RDS.
#
# Just the NMF step: PCA-elbow chooses num_factors (capped at --max_factors,
# floored at 2), then fastTopics::fit_poisson_nmf fits the model.
# No DE, no plots -- those are gene_level-only follow-up rules.
# No annotated/unannotated split here either -- that's the
# extract_unann_factors.R step.
#
# Input matrix.RDS schema (from ingest_snaptron / ingest_gene_matrix +
# classify_junctions):
#   list(counts        = junctions x samples integer matrix,
#        junction_info = tibble(junction_id, chrom, start, end, strand, annotated),
#        gene_name, gene_chrom, gene_strand)
#
# Output res.RDS:
#   - the fastTopics fit (rownames(F) = junction_id, rownames(L) = sample_id)
#   - OR a list(skipped = TRUE, reason = ..., gene_name = ...) pass-through
#     when upstream skipped or the PCA elbow yields no usable factor count.
#
# Replaces v0.1 elbow_cutoff (an inverted phrasing) with variance_explained:
#   v0.1: ideal_factors = length(cv) - sum(cv > (1 - elbow_cutoff))   # 0.01
#   v1.0: ideal_factors = length(cv) - sum(cv >  variance_explained)  # 0.99
# Math is identical when variance_explained = 1 - elbow_cutoff; this just
# names the knob after what it actually represents.

suppressPackageStartupMessages({
  library(optparse)
})
# Install fastTopics on first use (not packaged for conda-forge reliably).
if (!requireNamespace("fastTopics", quietly = TRUE)) {
  remotes::install_github("stephenslab/fastTopics",
                           upgrade = "never", dependencies = TRUE)
}
suppressPackageStartupMessages({
  library(fastTopics)
})

option_list <- list(
  make_option("--matrix", type = "character",
              help = "input matrix.RDS"),
  make_option("--output", type = "character",
              help = "output res.RDS"),
  make_option("--max_factors",        type = "integer", default = 10L,
              help = "upper bound on the per-gene k [%default]"),
  make_option("--variance_explained", type = "double",  default = 0.99,
              help = "cumulative PCA variance threshold for the elbow [%default]"),
  make_option("--threads", type = "integer", default = 1L,
              help = "fit_poisson_nmf control nc (parallel threads) [%default]"),
  make_option("--min_junctions_for_pca", type = "integer", default = 3L,
              help = "hard floor on rows/cols before running PCA [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$matrix), !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

mat <- readRDS(opt$matrix)

# Pass-through if upstream skipped.
if (isTRUE(mat$skipped)) {
  saveRDS(mat, opt$output)
  message(sprintf("  [pass-through %s] upstream skipped (%s)",
                   mat$gene_name, mat$reason))
  quit(save = "no", status = 0)
}

counts   <- mat$counts                       # junctions x samples (integer)
junc_ids <- mat$junction_info$junction_id
samples  <- colnames(counts)

# Sanity: post-canonicalise both dimensions should already be >= min_junctions_for_pca,
# but enforce again so an exotic matrix.RDS doesn't crash prcomp().
if (nrow(counts) < opt$min_junctions_for_pca ||
    ncol(counts) < opt$min_junctions_for_pca) {
  reason <- sprintf("matrix too small for PCA: %d junctions x %d samples (< %d)",
                    nrow(counts), ncol(counts), opt$min_junctions_for_pca)
  saveRDS(list(skipped = TRUE, reason = reason, gene_name = mat$gene_name),
           opt$output)
  message(sprintf("  [skip %s] %s", mat$gene_name, reason))
  quit(save = "no", status = 0)
}

# ---- PCA elbow chooses num_factors ---------------------------------------
# Normalise per sample so each sample's junction-count vector sums to 1
# ('equally weight' samples in PCA, same as v0.1). Then PCA over samples,
# with junctions as features.
sample_totals <- colSums(counts)
sample_totals[sample_totals == 0L] <- 1L      # defensive; shouldn't fire post-canonical
normalised    <- t(t(counts) / sample_totals) # junctions x samples (samples are columns, sum 1)
pc            <- prcomp(t(normalised))        # samples as observations
prop_var      <- as.numeric(summary(pc)$importance[2, ])
cumul_var     <- cumsum(prop_var)

ideal_factors <- length(prop_var) - sum(cumul_var > opt$variance_explained)
num_factors   <- max(2L, min(opt$max_factors, ideal_factors))

message(sprintf("  [pnmf %s] %d junctions x %d samples; PCA -> ideal_factors=%d, k=%d",
                 mat$gene_name, nrow(counts), ncol(counts),
                 ideal_factors, num_factors))

# ---- Poisson NMF ---------------------------------------------------------
dat <- t(as.matrix(counts))                   # samples x junctions, as fastTopics expects
res <- fit_poisson_nmf(dat, k = num_factors,
                        control = list(nc = opt$threads))

# Tag rows of F with the canonical junction_id (was "start-end" in v0.1).
rownames(res$F) <- junc_ids
rownames(res$L) <- samples

# Stash the gene_name on the res object so extract_unann_factors.R and
# downstream don't need to re-parse a path to recover it.
res$gene_name <- mat$gene_name

saveRDS(res, opt$output)
message(sprintf("  wrote %s", opt$output))
