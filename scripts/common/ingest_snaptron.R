#!/usr/bin/env Rscript
# Quetzal v1.0 -- snaptron TSV -> canonical per-gene matrix.RDS.
#
# Snaptron-format-specific parsing only; the post-parse pipeline (range/strand
# filter, normal-drop, leafcutter clustering, cluster/sample read floors,
# size checks) lives in scripts/common/_canonical.R so this script and
# ingest_gene_matrix.R stay in lockstep.
#
# Standalone invocation:
#   Rscript scripts/common/ingest_snaptron.R \
#       --input data/all_genes/chr5/snaptron_output/BRD9_snaptron.tsv \
#       --gene_name BRD9 --chr chr5 \
#       --gencode data/gencode_v46_granges.RDS \
#       --output output/gene_level/chr5/BRD9/matrix.RDS

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
              help = "snaptron TSV for one gene"),
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
              help = "junctions kept must fall within +/- N bp of the gene's gencode range [%default]"),
  make_option("--min_clust_read_count_avg", type = "double", default = 0,
              help = "clusters dropped if mean reads-per-sample < N [%default]"),
  make_option("--min_reads_per_sample_per_cluster", type = "double", default = 5,
              help = "samples dropped if mean reads-per-cluster < N [%default]"),

  # optional normal-sample exclusion (off by default in v1.0)
  make_option("--sample_metadata",       type = "character", default = NA_character_,
              help = "path to sample metadata TSV (needed for --exclude_normals)"),
  make_option("--exclude_normals",       type = "logical",   default = FALSE,
              action = "store_true",
              help = "drop samples whose --normal_filter_column matches --normal_filter_pattern"),
  make_option("--normal_filter_column",  type = "character", default = NA_character_,
              help = "metadata column to test (e.g. 'cgc_sample_sample_type' for TCGA)"),
  make_option("--normal_filter_pattern", type = "character", default = NA_character_,
              help = "regex pattern matching normal samples (e.g. 'Normal' for TCGA)"),
  make_option("--sample_id_column",      type = "character", default = "rail_id",
              help = "metadata column carrying sample IDs used as snaptron count keys [%default]"),

  # internal early-exit threshold (also enforced again in fit_pnmf.R)
  make_option("--min_junctions_for_pca", type = "integer", default = 3L,
              help = "hard floor on surviving junctions [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input), !is.null(opt$gene_name),
           !is.null(opt$chr),   !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

# ---- load shared helpers -------------------------------------------------

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

# ---- short-circuit: empty / tiny input -----------------------------------

if (!file.exists(opt$input) || file.size(opt$input) < 100L) {
  write_result(.skipped(opt$gene_name,
                          sprintf("input file missing or < 100 bytes (%s)", opt$input)),
                opt$output)
  quit(save = "no", status = 0)
}

# ---- snaptron TSV -> raw counts + junction_info --------------------------

gene_table <- read_tsv(opt$input, show_col_types = FALSE)
required   <- c("chromosome", "start", "end", "strand", "samples",
                 "samples_count", "annotated")
missing    <- setdiff(required, colnames(gene_table))
if (length(missing) > 0L) {
  write_result(.skipped(opt$gene_name,
                          sprintf("snaptron table missing required columns: %s",
                                  paste(missing, collapse = ", "))),
                opt$output)
  quit(save = "no", status = 0)
}

# Snaptron-specific cheap filter on the `samples_count` column. The
# stage-1 rowSums filter in _canonical.R catches the rest (and handles
# generic input that has no `samples_count`).
gene_table <- gene_table %>%
  filter(samples_count > opt$min_samples_per_junc) %>%
  distinct()

# Sample IDs: from metadata when supplied, else the union seen in the
# `samples` packed strings.
sample_meta <- NULL
if (!is.na(opt$sample_metadata) && nzchar(opt$sample_metadata)) {
  sample_meta <- read_tsv(opt$sample_metadata, show_col_types = FALSE)
  if (!opt$sample_id_column %in% colnames(sample_meta)) {
    write_result(.skipped(opt$gene_name,
                            sprintf("sample metadata missing --sample_id_column '%s'",
                                    opt$sample_id_column)),
                  opt$output)
    quit(save = "no", status = 0)
  }
  sample_ids <- as.character(sample_meta[[opt$sample_id_column]])
} else {
  sample_ids <- unique(unlist(lapply(gene_table$samples, function(s) {
    parts <- str_split(s, ",")[[1]]
    parts <- parts[nzchar(parts)]
    vapply(str_split(parts, ":"), `[`, character(1L), 1L)
  })))
}

counts <- matrix(0L,
                  nrow = nrow(gene_table),
                  ncol = length(sample_ids),
                  dimnames = list(NULL, sample_ids))

for (i in seq_len(nrow(gene_table))) {
  pre_split <- str_split(gene_table$samples[i], ",")[[1]]
  pre_split <- pre_split[nzchar(pre_split)]
  if (!length(pre_split)) next
  pairs     <- str_split_fixed(pre_split, ":", 2L)
  rail_ids  <- pairs[, 1]
  cs        <- as.integer(pairs[, 2])
  hit       <- rail_ids %in% colnames(counts)
  if (any(hit)) counts[i, rail_ids[hit]] <- cs[hit]
}

junction_info <- tibble(
  junction_id = paste0(gene_table$chromosome, ":",
                        gene_table$start, "-", gene_table$end, ":",
                        gene_table$strand),
  chrom       = gene_table$chromosome,
  start       = as.integer(gene_table$start),
  end         = as.integer(gene_table$end),
  strand      = gene_table$strand,
  annotated   = as.integer(gene_table$annotated)
)
rownames(counts) <- junction_info$junction_id

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
  message(sprintf("  [ok %s] %d junctions x %d samples",
                   opt$gene_name, nrow(result$counts), ncol(result$counts)))
  write_result(result, opt$output)
}
