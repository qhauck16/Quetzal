#!/usr/bin/env Rscript
# Quetzal v1.0 -- snaptron TSV -> canonical per-gene matrix.RDS.
#
# Input  : one snaptron-formatted TSV per gene. Required columns:
#            chromosome | start | end | strand | samples (string)
#            samples_count (int) | annotated (0/1)
#          The `samples` column is the snaptron-packed
#          `,railid1:count1,railid2:count2,...` string.
# Output : an RDS holding a list:
#            list(
#              counts        = matrix(int, junctions x samples),
#                                rownames = junction_id ("chr:start-end:strand"),
#                                colnames = sample_id
#              junction_info = tibble(junction_id, chrom, start, end, strand, annotated),
#              gene_name, gene_chrom, gene_strand
#            )
#          When the gene fails any filter, the RDS holds
#            list(skipped = TRUE, reason = "...")
#          so Snakemake still sees a file at the declared output path.
#
# Mirrors the data-prep half of v0.1 `tcga_LF_saving.R` (lines 73-241) but
# behind a clean CLI, with the v0.1 "drop TCGA normals" step gated behind an
# optional --exclude_normals / --normal_filter_* flag set.
#
# Standalone invocation:
#   Rscript ingest_snaptron.R \
#       --input data/all_genes/chr5/snaptron_output/BRD9_snaptron.tsv \
#       --gene_name BRD9 --chr chr5 \
#       --gencode data/hg38_granges.RDS \
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
              default = "data/hg38_granges.RDS",
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
  make_option("--sample_metadata",         type = "character", default = NA_character_,
              help = "path to sample metadata TSV (needed for --exclude_normals)"),
  make_option("--exclude_normals",         type = "logical",   default = FALSE,
              action = "store_true",
              help = "drop samples whose --normal_filter_column matches --normal_filter_pattern"),
  make_option("--normal_filter_column",    type = "character", default = NA_character_,
              help = "metadata column to test (e.g. 'cgc_sample_sample_type' for TCGA)"),
  make_option("--normal_filter_pattern",   type = "character", default = NA_character_,
              help = "regex pattern matching normal samples (e.g. 'Normal' for TCGA)"),
  make_option("--sample_id_column",        type = "character", default = "rail_id",
              help = "metadata column carrying the sample IDs used as snaptron count keys [%default]"),

  # internal early-exit threshold (also enforced again in fit_pnmf.R)
  make_option("--min_junctions_for_pca",   type = "integer",   default = 3L,
              help = "hard floor on surviving junctions; below this we skip the gene [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

stopifnot(!is.null(opt$input), !is.null(opt$gene_name),
           !is.null(opt$chr),   !is.null(opt$output))

dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

# ---- helpers --------------------------------------------------------------

skip <- function(reason) {
  message(sprintf("  [skip %s] %s", opt$gene_name, reason))
  saveRDS(list(skipped = TRUE, reason = reason, gene_name = opt$gene_name),
           opt$output)
  quit(save = "no", status = 0)
}

# ---- short-circuit: empty / tiny input ------------------------------------

if (!file.exists(opt$input) || file.size(opt$input) < 100L) {
  skip(sprintf("input file missing or < 100 bytes (%s)", opt$input))
}

# ---- gencode -> gene coords ----------------------------------------------

gencode    <- readRDS(opt$gencode)
gene_data  <- gencode[gencode$gene_name == opt$gene_name &
                        gencode$type     == "exon"]
if (length(gene_data) == 0L) {
  skip(sprintf("gene %s not found in gencode", opt$gene_name))
}

exon_starts <- start(gene_data)
exon_ends   <- end(gene_data)
lower_bound <- min(exon_starts)
upper_bound <- max(exon_ends)
gene_strand <- unique(as.character(strand(gene_data)))
gene_chrom  <- unique(as.character(seqnames(gene_data)))[1]

# ---- snaptron table -> junction filter -----------------------------------

gene_table <- read_tsv(opt$input, show_col_types = FALSE)
required   <- c("chromosome", "start", "end", "strand", "samples",
                 "samples_count", "annotated")
missing    <- setdiff(required, colnames(gene_table))
if (length(missing) > 0L) {
  skip(sprintf("snaptron table missing required columns: %s",
               paste(missing, collapse = ", ")))
}

gene_table <- gene_table %>%
  filter(samples_count > opt$min_samples_per_junc) %>%
  distinct() %>%
  filter(start > (lower_bound - opt$gene_range_bound)) %>%
  filter(end   < (upper_bound + opt$gene_range_bound)) %>%
  filter(strand %in% gene_strand)

if (nrow(gene_table) < opt$min_junctions_for_pca) {
  skip(sprintf("only %d junctions after sample/strand/range filter (< %d)",
               nrow(gene_table), opt$min_junctions_for_pca))
}

# ---- parse snaptron `samples` column -> junction x sample count matrix ---

# Sample IDs come from sample metadata when supplied, otherwise harvested
# directly from snaptron's samples strings (union across all junctions).
if (!is.na(opt$sample_metadata) && nzchar(opt$sample_metadata)) {
  sample_meta <- read_tsv(opt$sample_metadata, show_col_types = FALSE)
  if (!opt$sample_id_column %in% colnames(sample_meta)) {
    skip(sprintf("sample metadata missing --sample_id_column '%s'",
                 opt$sample_id_column))
  }
  sample_ids <- as.character(sample_meta[[opt$sample_id_column]])
} else {
  sample_meta <- NULL
  sample_ids  <- unique(unlist(lapply(gene_table$samples, function(s) {
    parts <- str_split(s, ",")[[1]]
    parts <- parts[nzchar(parts)]
    vapply(str_split(parts, ":"), `[`, character(1L), 1L)
  })))
}

just_counts <- matrix(0L,
                       nrow = nrow(gene_table),
                       ncol = length(sample_ids),
                       dimnames = list(NULL, sample_ids))

for (i in seq_len(nrow(gene_table))) {
  pre_split <- str_split(gene_table$samples[i], ",")[[1]]
  pre_split <- pre_split[nzchar(pre_split)]
  if (!length(pre_split)) next
  pairs     <- str_split_fixed(pre_split, ":", 2L)
  rail_ids  <- pairs[, 1]
  counts    <- as.integer(pairs[, 2])
  hit       <- rail_ids %in% colnames(just_counts)
  if (any(hit)) {
    just_counts[i, rail_ids[hit]] <- counts[hit]
  }
}

# ---- optional normal-sample filter ----------------------------------------

if (isTRUE(opt$exclude_normals)) {
  if (is.null(sample_meta))               skip("--exclude_normals set but no --sample_metadata supplied")
  if (is.na(opt$normal_filter_column))    skip("--exclude_normals set but no --normal_filter_column supplied")
  if (is.na(opt$normal_filter_pattern))   skip("--exclude_normals set but no --normal_filter_pattern supplied")
  if (!opt$normal_filter_column %in% colnames(sample_meta)) {
    skip(sprintf("sample metadata missing --normal_filter_column '%s'",
                 opt$normal_filter_column))
  }
  normal_mask <- grepl(opt$normal_filter_pattern,
                        sample_meta[[opt$normal_filter_column]])
  normal_ids  <- as.character(sample_meta[[opt$sample_id_column]][normal_mask])
  just_counts <- just_counts[, !colnames(just_counts) %in% normal_ids, drop = FALSE]
}

# ---- re-apply min_samples_per_junc post sample-set change ----------------

keep <- rowSums(just_counts > 0L) >= opt$min_samples_per_junc
gene_table  <- gene_table[keep, , drop = FALSE]
just_counts <- just_counts[keep, , drop = FALSE]

if (nrow(gene_table) < opt$min_junctions_for_pca) {
  skip(sprintf("only %d junctions after sample-set rebuild (< %d)",
               nrow(gene_table), opt$min_junctions_for_pca))
}

# ---- leafcutter-esque clustering -----------------------------------------
# Group junctions sharing a `start`, then merge clusters that share an `end`
# with an earlier cluster. Mirrors v0.1 tcga_LF_saving.R lines 177-192.

clustering <- gene_table %>% mutate(.row_idx = row_number()) %>%
                              group_by(start)    %>%
                              mutate(cluster_id = cur_group_id()) %>%
                              ungroup()
clusters <- clustering$cluster_id

if (length(clusters) > 1L) {
  for (i in 2:max(clusters)) {
    ends_i <- clustering$end[clusters == i]
    if (any(ends_i %in% clustering$end[clusters < i])) {
      match_idx <- which(clustering$end %in% ends_i & clusters < i)[1]
      clusters[clusters == i] <- clusters[match_idx]
    }
  }
}
gene_table$cluster <- clusters

# ---- per-cluster sums + cluster + sample filtering -----------------------

clust_sums <- rowsum(just_counts, group = clusters, reorder = FALSE)

keep_clust   <- (rowMeans(clust_sums) >= opt$min_clust_read_count_avg)
keep_cluster_ids <- as.integer(rownames(clust_sums)[keep_clust])
if (!any(keep_clust)) {
  skip(sprintf("no clusters passed min_clust_read_count_avg = %s",
               opt$min_clust_read_count_avg))
}

j_keep <- gene_table$cluster %in% keep_cluster_ids
gene_table  <- gene_table[j_keep, , drop = FALSE]
just_counts <- just_counts[j_keep, , drop = FALSE]
clust_sums  <- clust_sums[keep_clust, , drop = FALSE]

keep_sample <- (colMeans(clust_sums) >= opt$min_reads_per_sample_per_cluster)
if (!any(keep_sample)) {
  skip(sprintf("no samples passed min_reads_per_sample_per_cluster = %s",
               opt$min_reads_per_sample_per_cluster))
}
just_counts <- just_counts[, keep_sample, drop = FALSE]

# ---- final dead-row sweep + size check -----------------------------------

keep_final  <- rowSums(just_counts) > 0L
just_counts <- just_counts[keep_final, , drop = FALSE]
gene_table  <- gene_table[keep_final, , drop = FALSE]

if (nrow(just_counts) < opt$min_junctions_for_pca ||
    ncol(just_counts) < opt$min_junctions_for_pca) {
  skip(sprintf("final matrix %d junctions x %d samples (< %d)",
               nrow(just_counts), ncol(just_counts), opt$min_junctions_for_pca))
}

# ---- canonical junction info ---------------------------------------------

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
rownames(just_counts) <- junction_info$junction_id

out <- list(
  counts        = just_counts,
  junction_info = junction_info,
  gene_name     = opt$gene_name,
  gene_chrom    = gene_chrom,
  gene_strand   = paste(gene_strand, collapse = ",")
)

message(sprintf("  [ok %s] %d junctions x %d samples written to %s",
                 opt$gene_name, nrow(just_counts), ncol(just_counts), opt$output))
saveRDS(out, opt$output)
