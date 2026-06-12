#!/usr/bin/env Rscript
# Quetzal v1.0 -- fill in matrix.RDS$junction_info$annotated by comparing
# each junction's (start, end) against the gencode-derived annotated intron
# set for the gene.
#
# This step is only needed for `input_format: gene_matrix`. The snaptron
# ingest path takes `annotated` straight from snaptron's column, so the
# gene_level / genome_wide Snakefiles skip this rule entirely for that
# input format.
#
# Definition: a junction (start, end) is annotated if some gencode
# transcript of this gene has consecutive exons whose intron coordinates
# (upstream exon end + 1, downstream exon start - 1) match (start, end)
# exactly. Coords are 1-based inclusive on both ends.

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tibble)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--matrix",  type = "character",
              help = "input matrix.RDS (junction_info$annotated may be NA)"),
  make_option("--gencode", type = "character",
              default = "data/hg38_granges.RDS",
              help = "gencode source (.RDS) [%default]"),
  make_option("--output",  type = "character",
              help = "output matrix.RDS with junction_info$annotated populated")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$matrix), !is.null(opt$output))
dir.create(dirname(opt$output), showWarnings = FALSE, recursive = TRUE)

mat <- readRDS(opt$matrix)

# Pass through skipped stubs unchanged.
if (isTRUE(mat$skipped)) {
  saveRDS(mat, opt$output)
  message(sprintf("  [pass-through %s] upstream skipped (%s)",
                   mat$gene_name, mat$reason))
  quit(save = "no", status = 0)
}

gencode   <- readRDS(opt$gencode)
gene_data <- gencode[gencode$gene_name == mat$gene_name &
                        gencode$type     == "exon"]

if (length(gene_data) == 0L) {
  mat$junction_info$annotated <- 0L
  saveRDS(mat, opt$output)
  message(sprintf("  [annot %s] no gencode exons; all %d junctions marked unannotated",
                   mat$gene_name, nrow(mat$junction_info)))
  quit(save = "no", status = 0)
}

# Per-transcript: sort exons by genomic start, derive intron coords from
# consecutive gaps. Genomic ordering works on both strands because
# junctions are genomic-coord too.
exon_df <- as.data.frame(gene_data) %>%
  dplyr::select(transcript_id, start, end) %>%
  arrange(transcript_id, start)

annotated_introns <- exon_df %>%
  group_by(transcript_id) %>%
  mutate(intron_start = end + 1L,
          intron_end   = lead(start) - 1L) %>%
  filter(!is.na(intron_end), intron_start <= intron_end) %>%
  ungroup() %>%
  distinct(intron_start, intron_end) %>%
  mutate(key = paste0(intron_start, "-", intron_end)) %>%
  pull(key)

mat$junction_info <- mat$junction_info %>%
  mutate(annotated = as.integer(
    paste0(start, "-", end) %in% annotated_introns
  ))

n_ann <- sum(mat$junction_info$annotated == 1L)
n_tot <- nrow(mat$junction_info)
message(sprintf("  [annot %s] %d / %d junctions annotated",
                 mat$gene_name, n_ann, n_tot))
saveRDS(mat, opt$output)
