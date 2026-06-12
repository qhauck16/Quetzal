#!/usr/bin/env Rscript
# Quetzal v1.0 helper -- leafcutter-style splice-junction count matrix
# -> per-gene wide-TSV files in the format ingest_gene_matrix.R consumes.
#
# Input file (CSV or TSV; auto-detected from extension):
#   col 1 ("title"): leafcutter junction ID "chr:start:end:clu_<id>_<strand>"
#                    (start is 0-based / BED-style; converted to 1-based here)
#   col 2..N       : sample IDs (header row), integer read counts in each cell.
#
# Per-gene filter (mirrors v0.1 tcga_LF_saving.R):
#   1. gencode-derived gene range [lower, upper] from `type == "exon"` rows
#      whose `gene_name` matches.
#   2. keep junctions where
#        chrom == gene chrom
#        strand == gene strand
#        start > (lower - gene_range_bound)
#        end   < (upper + gene_range_bound)
#
# Output (per gene that has >=1 surviving junction):
#   <output_dir>/<chr>/<gene>.tsv
#
# with the header row:
#   junction_id  S1  S2  S3  ...
# and junction_id formatted as "chr:start-end:strand" (1-based, matching
# the canonical convention ingest_gene_matrix.R expects).
#
# CLI:
#   Rscript scripts/helpers/leafcutter_to_gene_matrices.R \
#       --input  /path/to/leafcutter_counts.csv \
#       --gencode data/hg38_granges.RDS \
#       --output_dir data/all_genes \
#       --gene_range_bound 2000

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--input",            type = "character",
              help = "leafcutter-style junction count matrix (CSV or TSV)"),
  make_option("--gencode",          type = "character",
              default = "data/hg38_granges.RDS",
              help = "gencode source (.RDS) [%default]"),
  make_option("--output_dir",       type = "character",
              default = "data/all_genes",
              help = "root output dir; per-chr subdirs created here [%default]"),
  make_option("--gene_range_bound", type = "integer", default = 2000L,
              help = "bp padding around gene's exon range [%default]"),
  make_option("--min_junctions",    type = "integer", default = 1L,
              help = "only write per-gene TSV if it has at least N junctions [%default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$input))
dir.create(opt$output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- read input ----------------------------------------------------------

sep <- if (grepl("\\.csv$", opt$input, ignore.case = TRUE)) "," else "\t"
message(sprintf("reading %s (sep='%s') ...", opt$input, sep))
dt <- fread(opt$input, sep = sep, header = TRUE, showProgress = TRUE)
message(sprintf("  loaded %s junctions x %s samples",
                 format(nrow(dt), big.mark=","),
                 format(ncol(dt) - 1L, big.mark=",")))

id_col   <- colnames(dt)[1]
samples  <- colnames(dt)[-1L]
ids_raw  <- dt[[id_col]]

# ---- parse leafcutter junction IDs --------------------------------------
# Format: chr:start:end:clu_<id>_<strand>
# Start is BED-style 0-based; convert to 1-based by adding 1.

parts <- str_split_fixed(ids_raw, ":", 4L)
chrom         <- parts[, 1]
start_1based  <- as.integer(parts[, 2]) + 1L
end_coord     <- as.integer(parts[, 3])
cluster_field <- parts[, 4]
strand        <- str_extract(cluster_field, "[+\\-]$")

bad <- is.na(start_1based) | is.na(end_coord) | is.na(strand) |
        !nzchar(chrom)
if (any(bad)) {
  message(sprintf("  WARNING: %d junction IDs failed to parse and will be dropped (e.g. '%s')",
                   sum(bad), ids_raw[which(bad)[1]]))
  dt           <- dt[!bad]
  chrom        <- chrom[!bad]
  start_1based <- start_1based[!bad]
  end_coord    <- end_coord[!bad]
  strand       <- strand[!bad]
  ids_raw      <- ids_raw[!bad]
}

# canonical Quetzal v1.0 junction_id = "chr:start-end:strand"
junction_id <- sprintf("%s:%d-%d:%s", chrom, start_1based, end_coord, strand)

# ---- gencode -> gene_name -> (chr, strand, lower, upper) ----------------

message(sprintf("loading gencode from %s ...", opt$gencode))
gencode  <- readRDS(opt$gencode)
exons    <- gencode[gencode$type == "exon"]
gene_df  <- as.data.frame(exons) %>%
  group_by(gene_name) %>%
  summarise(
    chr    = as.character(seqnames)[1],
    strand = as.character(strand)[1],
    lower  = min(start),
    upper  = max(end),
    .groups = "drop"
  ) %>%
  filter(!is.na(gene_name), nzchar(gene_name))
message(sprintf("  %d unique gene_names", nrow(gene_df)))

# ---- per-chr index over the junction set --------------------------------

idx_by_chr <- split(seq_along(chrom), chrom)

# ---- per-gene slice + write ---------------------------------------------

message(sprintf("writing per-gene TSVs under %s/ ...", opt$output_dir))
n_written <- 0L
n_skipped <- 0L

for (i in seq_len(nrow(gene_df))) {
  gn       <- gene_df$gene_name[i]
  g_chr    <- gene_df$chr[i]
  g_strand <- gene_df$strand[i]
  g_lower  <- gene_df$lower[i]  - opt$gene_range_bound
  g_upper  <- gene_df$upper[i]  + opt$gene_range_bound

  idx <- idx_by_chr[[g_chr]]
  if (is.null(idx) || !length(idx)) { n_skipped <- n_skipped + 1L; next }

  keep <- idx[
    strand[idx]       == g_strand &
    start_1based[idx] >  g_lower  &
    end_coord[idx]    <  g_upper
  ]
  if (length(keep) < opt$min_junctions) {
    n_skipped <- n_skipped + 1L
    next
  }

  out_dir <- file.path(opt$output_dir, g_chr)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_path <- file.path(out_dir, paste0(gn, ".tsv"))

  # build wide tibble: junction_id + sample columns
  sub <- dt[keep]
  sub[[id_col]] <- junction_id[keep]
  setnames(sub, id_col, "junction_id")

  fwrite(sub, out_path, sep = "\t")
  n_written <- n_written + 1L

  if (n_written %% 500L == 0L) {
    message(sprintf("  ... %d / %d genes written",
                     n_written, nrow(gene_df)))
  }
}

message(sprintf("done: wrote %d per-gene TSVs; skipped %d (no surviving junctions / unknown chrom)",
                 n_written, n_skipped))
