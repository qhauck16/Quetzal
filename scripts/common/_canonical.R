# Shared post-parse pipeline for the per-gene canonical matrix.RDS.
# Used by ingest_snaptron.R and ingest_gene_matrix.R; not directly executable.
#
# Expected callers' raw inputs:
#   counts        : integer matrix, junctions x samples
#                   rownames = junction_id ("chr:start-end:strand")
#                   colnames = sample_id (string)
#   junction_info : tibble(junction_id, chrom, start, end, strand, annotated)
#                   `annotated` may be NA_integer_ when caller can't fill it
#                   (gene_matrix path), in which case classify_junctions.R
#                   fills it later.
#   gene_name     : character(1)
#   gencode       : the readRDS()'d gencode GRanges
#   opt           : list with these knobs (parsed from each caller's CLI):
#                     min_samples_per_junc, gene_range_bound,
#                     min_clust_read_count_avg,
#                     min_reads_per_sample_per_cluster,
#                     min_junctions_for_pca,
#                     sample_metadata (data.frame or NULL),
#                     exclude_normals (logical),
#                     normal_filter_column, normal_filter_pattern,
#                     sample_id_column
#
# Returns either:
#   list(skipped = TRUE, reason = "...", gene_name = ...)   on early-exit
# or:
#   list(counts, junction_info, gene_name, gene_chrom, gene_strand)
#
# Callers saveRDS() the result themselves.

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(GenomicRanges)
})

.skipped <- function(gene_name, reason) {
  message(sprintf("  [skip %s] %s", gene_name, reason))
  list(skipped = TRUE, reason = reason, gene_name = gene_name)
}

# Returns list(lower, upper, strand, chrom) or NULL if gene absent from gencode.
.resolve_gene_coords <- function(gencode, gene_name) {
  gene_data <- gencode[gencode$gene_name == gene_name &
                          gencode$type     == "exon"]
  if (length(gene_data) == 0L) return(NULL)
  list(
    lower  = min(start(gene_data)),
    upper  = max(end(gene_data)),
    strand = unique(as.character(strand(gene_data))),
    chrom  = unique(as.character(seqnames(gene_data)))[1]
  )
}

# Leafcutter-esque clustering: group junctions sharing a `start`, then merge
# clusters that share an `end` with an earlier cluster. Mirrors v0.1
# tcga_LF_saving.R lines 177-192.
.cluster_junctions <- function(junction_info) {
  starts <- junction_info$start
  ends   <- junction_info$end
  ord    <- order(starts)
  ids    <- as.integer(factor(starts[ord], levels = unique(starts[ord])))
  # ids[i] is the cluster id for junction `ord[i]`. Walk in cluster order and
  # merge if an end matches an earlier cluster's end.
  if (length(ids) > 1L) {
    for (i in 2:max(ids)) {
      ends_i <- ends[ord][ids == i]
      if (any(ends_i %in% ends[ord][ids < i])) {
        match_idx <- which(ends[ord] %in% ends_i & ids < i)[1]
        ids[ids == i] <- ids[match_idx]
      }
    }
  }
  # invert the permutation
  out <- integer(length(ids))
  out[ord] <- ids
  out
}

canonicalise <- function(counts, junction_info, gene_name, gencode, opt) {

  stopifnot(is.matrix(counts), is.data.frame(junction_info))
  stopifnot(nrow(counts) == nrow(junction_info))
  if (is.null(rownames(counts))) rownames(counts) <- junction_info$junction_id

  # ---- resolve gene coords + range/strand junction filter ---------------

  gene_coords <- .resolve_gene_coords(gencode, gene_name)
  if (is.null(gene_coords)) {
    return(.skipped(gene_name, "gene not found in gencode"))
  }

  rb <- opt$gene_range_bound
  keep <- (junction_info$start > (gene_coords$lower - rb)) &
          (junction_info$end   < (gene_coords$upper + rb)) &
          (junction_info$strand %in% gene_coords$strand)
  junction_info <- junction_info[keep, , drop = FALSE]
  counts        <- counts[keep, , drop = FALSE]

  # ---- min_samples_per_junc (stage 1) -----------------------------------

  keep <- rowSums(counts > 0L) >= opt$min_samples_per_junc
  junction_info <- junction_info[keep, , drop = FALSE]
  counts        <- counts[keep, , drop = FALSE]

  if (nrow(junction_info) < opt$min_junctions_for_pca) {
    return(.skipped(gene_name,
                     sprintf("%d junctions after stage-1 filter (< %d)",
                             nrow(junction_info), opt$min_junctions_for_pca)))
  }

  # ---- optional normal-sample filter ------------------------------------

  if (isTRUE(opt$exclude_normals)) {
    if (is.null(opt$sample_metadata))            return(.skipped(gene_name, "--exclude_normals set but no sample_metadata"))
    if (is.na(opt$normal_filter_column))         return(.skipped(gene_name, "--exclude_normals set but no normal_filter_column"))
    if (is.na(opt$normal_filter_pattern))        return(.skipped(gene_name, "--exclude_normals set but no normal_filter_pattern"))
    if (!opt$normal_filter_column %in% colnames(opt$sample_metadata)) {
      return(.skipped(gene_name,
                       sprintf("metadata missing column '%s'",
                               opt$normal_filter_column)))
    }
    nm <- grepl(opt$normal_filter_pattern,
                 opt$sample_metadata[[opt$normal_filter_column]])
    nm_ids <- as.character(opt$sample_metadata[[opt$sample_id_column]][nm])
    counts <- counts[, !colnames(counts) %in% nm_ids, drop = FALSE]
  }

  # ---- min_samples_per_junc (stage 2, post-sample-set change) ----------

  keep <- rowSums(counts > 0L) >= opt$min_samples_per_junc
  junction_info <- junction_info[keep, , drop = FALSE]
  counts        <- counts[keep, , drop = FALSE]

  if (nrow(junction_info) < opt$min_junctions_for_pca) {
    return(.skipped(gene_name,
                     sprintf("%d junctions after stage-2 filter (< %d)",
                             nrow(junction_info), opt$min_junctions_for_pca)))
  }

  # ---- leafcutter clustering + cluster/sample read floors --------------

  clusters    <- .cluster_junctions(junction_info)
  clust_sums  <- rowsum(counts, group = clusters, reorder = FALSE)
  keep_clust  <- (rowMeans(clust_sums) >= opt$min_clust_read_count_avg)
  if (!any(keep_clust)) {
    return(.skipped(gene_name,
                     sprintf("no clusters passed min_clust_read_count_avg = %s",
                             opt$min_clust_read_count_avg)))
  }
  keep_cluster_ids <- as.integer(rownames(clust_sums)[keep_clust])

  j_keep        <- clusters %in% keep_cluster_ids
  junction_info <- junction_info[j_keep, , drop = FALSE]
  counts        <- counts[j_keep, , drop = FALSE]
  clust_sums    <- clust_sums[keep_clust, , drop = FALSE]

  keep_sample <- (colMeans(clust_sums) >= opt$min_reads_per_sample_per_cluster)
  if (!any(keep_sample)) {
    return(.skipped(gene_name,
                     sprintf("no samples passed min_reads_per_sample_per_cluster = %s",
                             opt$min_reads_per_sample_per_cluster)))
  }
  counts <- counts[, keep_sample, drop = FALSE]

  # ---- final dead-row sweep + size check -------------------------------

  keep <- rowSums(counts) > 0L
  counts        <- counts[keep, , drop = FALSE]
  junction_info <- junction_info[keep, , drop = FALSE]

  if (nrow(counts) < opt$min_junctions_for_pca ||
      ncol(counts) < opt$min_junctions_for_pca) {
    return(.skipped(gene_name,
                     sprintf("final matrix %d junctions x %d samples (< %d)",
                             nrow(counts), ncol(counts), opt$min_junctions_for_pca)))
  }

  list(
    counts        = counts,
    junction_info = junction_info,
    gene_name     = gene_name,
    gene_chrom    = gene_coords$chrom,
    gene_strand   = paste(gene_coords$strand, collapse = ",")
  )
}

# Caller boilerplate: resolves the helper's own directory so each ingest
# script can `source(file.path(.canonical_dir, "_canonical.R"))` reliably
# regardless of where snakemake invokes it from.
.find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa   <- regmatches(args, regexpr("(?<=^--file=).+", args, perl = TRUE))
  if (length(fa)) return(dirname(normalizePath(fa[1])))
  stop("could not resolve script directory (run via Rscript)")
}
