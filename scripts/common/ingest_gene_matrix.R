#!/usr/bin/env Rscript
# Quetzal v1.0 -- user-supplied wide-TSV per-gene junction count matrix
# -> canonical per-gene matrix.RDS.
#
# Expected input format (one TSV per gene):
#   first column          : junction_id, "chr:start-end:strand"
#                            (header may be empty or "junction_id")
#   remaining columns     : sample IDs (header row)
#   cells                 : integer read counts
#
# The `annotated` field of junction_info is left NA at this stage --
# classify_junctions.R fills it from the gencode model in a follow-up
# Snakemake rule.

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
  library(GenomicRanges)
})

# ---- CLI ------------------------------------------------------------------

option_list <- list(
  make_option("--input",     type = "character",
              help = "wide TSV for one gene (junctions x samples)"),
  make_option("--gene_name", type = "character",
              help = "gene name (matches gencode gene_name)"),
  make_option("--chr",       type = "character",
              help = "chromosome (used for output path + sanity)"),
  make_option("--gencode",   type = "character",
              default = "data/gencode_v46_granges.RDS",
              help = "gencode source (.RDS) [%default]"),
  make_option("--output",    type = "character",
              help = "matrix.RDS output path"),

  # filter knobs
  make_option("--min_samples_per_junc", type = "integer", default = 10L,
              help = "junctions must be supported in >= N samples [%default]"),
  make_option("--gene_range_bound",     type = "integer", default = 2000L,
              help = "junctions within +/- N bp of the gene's gencode range [%default]"),
  make_option("--min_clust_read_count_avg", type = "double", default = 0,
              help = "clusters dropped if mean reads-per-sample < N [%default]"),
  make_option("--min_reads_per_sample_per_cluster", type = "double", default = 5,
              help = "samples dropped if mean reads-per-cluster < N [%default]"),

  # optional normal-sample exclusion (off by default in v1.0)
  make_option("--sample_metadata",       type = "character", default = NA_character_,
              help = "path to sample metadata TSV"),
  make_option("--exclude_normals",       type = "logical",   default = FALSE,
              action = "store_true",
              help = "drop samples whose --normal_filter_column matches --normal_filter_pattern"),
  make_option("--normal_filter_column",  type = "character", default = NA_character_,
              help = "metadata column to test"),
  make_option("--normal_filter_pattern", type = "character", default = NA_character_,
              help = "regex pattern matching normal samples"),
  make_option("--sample_id_column",      type = "character", default = "sample_id",
              help = "metadata column carrying sample IDs (matches TSV column headers) [%default]"),

  make_option("--min_junctions_for_pca", type = "integer", default = 3L,
              help = "hard floor on surviving junctions [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$gene_name),
           !is.null(opt$chr),   !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

# ---- shared helpers -------------------------------------------------------

.find_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  fa   <- regmatches(args, regexpr("(?<=^--file=).+", args, perl = TRUE))
  if (length(fa)) return(dirname(normalizePath(fa[1])))
  stop("could not resolve script directory (run via Rscript)")
}
source(file.path(.find_script_dir(), "_canonical.R"))

write_result <- function(result, output) {
  saveRDS(result, output)
  message(sprintf("  wrote %s", output))
}

# ---- short-circuit -------------------------------------------------------

if (!file.exists(opt$input) || file.size(opt$input) < 50L) {
  write_result(.skipped(opt$gene_name,
                          sprintf("input file missing or < 50 bytes (%s)", opt$input)),
                opt$output)
  quit(save = "no", status = 0)
}

# ---- parse wide TSV -> raw counts + junction_info ------------------------

tbl <- read_tsv(opt$input, show_col_types = FALSE)
if (ncol(tbl) < 2L) {
  write_result(.skipped(opt$gene_name,
                          "input TSV has < 2 columns (need junction_id + at least one sample)"),
                opt$output)
  quit(save = "no", status = 0)
}

junction_id <- as.character(tbl[[1]])
sample_ids  <- colnames(tbl)[-1L]

counts <- as.matrix(tbl[, -1L, drop = FALSE])
mode(counts) <- "integer"
rownames(counts) <- junction_id
colnames(counts) <- sample_ids

# Parse junction_id "chr:start-end:strand" into components.
m <- regmatches(junction_id,
                 regexec("^([^:]+):([0-9]+)-([0-9]+):([+\\-\\*])$", junction_id))
bad <- vapply(m, function(x) length(x) == 0L, logical(1L))
if (any(bad)) {
  write_result(.skipped(opt$gene_name,
                          sprintf("%d junction IDs don't match 'chr:start-end:strand' (e.g. '%s')",
                                  sum(bad), junction_id[which(bad)[1]])),
                opt$output)
  quit(save = "no", status = 0)
}
parts <- do.call(rbind, m)
junction_info <- tibble(
  junction_id = junction_id,
  chrom       = parts[, 2],
  start       = as.integer(parts[, 3]),
  end         = as.integer(parts[, 4]),
  strand      = parts[, 5],
  annotated   = NA_integer_       # filled in by classify_junctions.R
)

# Optional sample metadata load (only used by canonicalise() for the
# normal filter; gene_matrix input doesn't need it otherwise).
sample_meta <- NULL
if (!is.na(opt$sample_metadata) && nzchar(opt$sample_metadata)) {
  sample_meta <- read_tsv(opt$sample_metadata, show_col_types = FALSE)
  if (!opt$sample_id_column %in% colnames(sample_meta)) {
    write_result(.skipped(opt$gene_name,
                            sprintf("metadata missing --sample_id_column '%s'",
                                    opt$sample_id_column)),
                  opt$output)
    quit(save = "no", status = 0)
  }
}

# ---- delegate to shared post-parse pipeline ------------------------------

gencode <- readRDS(opt$gencode)
result  <- canonicalise(
  counts        = counts,
  junction_info = junction_info,
  gene_name     = opt$gene_name,
  gencode       = gencode,
  opt = list(
    min_samples_per_junc             = opt$min_samples_per_junc,
    gene_range_bound                 = opt$gene_range_bound,
    min_clust_read_count_avg         = opt$min_clust_read_count_avg,
    min_reads_per_sample_per_cluster = opt$min_reads_per_sample_per_cluster,
    min_junctions_for_pca            = opt$min_junctions_for_pca,
    sample_metadata                  = sample_meta,
    exclude_normals                  = opt$exclude_normals,
    normal_filter_column             = opt$normal_filter_column,
    normal_filter_pattern            = opt$normal_filter_pattern,
    sample_id_column                 = opt$sample_id_column
  )
)

if (isTRUE(result$skipped)) {
  write_result(result, opt$output)
} else {
  message(sprintf("  [ok %s] %d junctions x %d samples (annotated NA, pending classify_junctions.R)",
                   opt$gene_name, nrow(result$counts), ncol(result$counts)))
  write_result(result, opt$output)
}
