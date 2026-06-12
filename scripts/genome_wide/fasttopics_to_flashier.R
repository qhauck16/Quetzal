#!/usr/bin/env Rscript
# Quetzal v0.1 - fastTopics -> softImpute -> flashier end-to-end
#
# Takes the gene-level fastTopics output directories produced by the
# genome_wide leafcutter snakefile (lf_Snakefile) and produces the final
# genome-wide multinomial flashier object used by all downstream
# (factor characterisation, NMD investigation, diff-splicing) analyses.
#
# Steps:
#   1. Iterate every (chr, gene) FastTopics_output/res.RDS:
#        - keep genes broadly expressed across TCGA (>= 80% of samples)
#        - keep local factors that load primarily on UNANNOTATED junctions
#          (max-loading on un-annotated > 0.5 * max-loading overall)
#        - convert to multinomial L matrix and rename columns
#          `<gene>.<kN>`.
#   2. full-join those per-gene L matrices into one big sample x feature
#      matrix and write it as `multinom_tcga_just_unn_just_factors_of_interest.tsv`.
#   3. Load TCGA RIN scores and per-chr productive/unproductive read
#      ratios; flag samples failing any of six QC tests (z > 3 / z < -3
#      on appropriate metrics, RIN < 5, junction_count < 125k).
#   4. Drop normal + QC-failing samples, fill missing values with
#      softImpute, scale, then fit a 300-factor `flashier::flash()` with
#      point-exponential + point-laplace EBNM priors.
#   5. Write the fitted `flash` object to disk.
#
# Usage (example):
#   Rscript scripts/genome_wide/fasttopics_to_flashier.R \
#       --gene_dir   ../260116_filters \
#       --snaptron_root ../../all_genes \
#       --tcga_meta  data/tcga_v2_samples.tsv \
#       --analyte    data/analyte.tsv \
#       --pu_dir     data/productive_unproductive \
#       --features_out softimpute_features.tsv \
#       --flash_out  softimpute_flash_300_qc_filtered.RDS
#
# All paths default to the layout shipped in the Quetzal v0.1 data/
# directory; override any of them with the matching CLI flag.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(tidyverse)
  library(fastTopics)
  library(flashier)
  library(softImpute)
})

option_list <- list(
  make_option("--gene_dir",       default = "../260116_filters",
              help = "root directory holding chr<N>/FastTopics_output/<GENE>/res.RDS"),
  make_option("--snaptron_root",  default = "../../all_genes",
              help = "root directory holding chr<N>/snaptron_output/<GENE>_snaptron.tsv"),
  make_option("--tcga_meta",      default = "data/tcga_v2_samples.tsv",
              help = "TCGA v2 sample metadata TSV"),
  make_option("--analyte",        default = "data/analyte.tsv",
              help = "TCGA analyte table (carries RIN scores)"),
  make_option("--pu_dir",         default = "data/productive_unproductive",
              help = "directory of chr<N>_productive_unproductive.tsv files"),
  make_option("--sample_fraction", type = "numeric", default = 0.8,
              help = "min fraction of TCGA samples a gene must cover [default %default]"),
  make_option("--features_out",   default = "softimpute_features.tsv",
              help = "intermediate per-sample factor-of-interest matrix"),
  make_option("--flash_out",      default = "softimpute_flash_300_qc_filtered.RDS",
              help = "final flashier object"),
  make_option("--greedy_kmax",    type = "integer", default = 300L,
              help = "flashier greedy_Kmax [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- 1. Per-gene unannotated-factor extraction -----------------------------

tcga_metadata    <- fread(opt$tcga_meta, data.table = FALSE)
normal_rail_ids  <- tcga_metadata$rail_id[
  grepl("Normal", tcga_metadata$cgc_sample_sample_type)]

unfilt <- data.frame(rail_id = as.character(
  tcga_metadata$rail_id[!tcga_metadata$rail_id %in% normal_rail_ids]),
  stringsAsFactors = FALSE)
sample_threshold <- opt$sample_fraction * nrow(unfilt)

chrs <- c(as.character(1:22), "X")
for (chr in chrs) {
  chr_dir   <- file.path(opt$gene_dir, paste0("chr", chr), "FastTopics_output")
  snap_dir  <- file.path(opt$snaptron_root, paste0("chr", chr), "snaptron_output")
  if (!dir.exists(chr_dir)) next
  genes <- basename(list.dirs(chr_dir, recursive = FALSE))

  per_gene_L <- list()
  unfilt_input <- data.frame(rail_id = as.character(tcga_metadata$rail_id),
                              stringsAsFactors = FALSE)

  k <- 0L
  for (gene in genes) {
    k <- k + 1L
    if (k %% 100L == 0L) message(sprintf("  chr%s  %d/%d", chr, k, length(genes)))

    res_path  <- file.path(chr_dir,  gene, "res.RDS")
    snap_path <- file.path(snap_dir, paste0(gene, "_snaptron.tsv"))
    if (file.size(res_path) <= 100) next

    res <- readRDS(res_path)
    if (nrow(res$L) < sample_threshold) next
    snap <- fread(snap_path) %>% filter(annotated == 1)
    annotated_juncs <- paste0(snap$start, "-", snap$end)

    unn_F <- as.matrix(res$F[!rownames(res$F) %in% annotated_juncs, , drop = FALSE])
    if (nrow(unn_F) == 0L) next

    factor_max_unn   <- apply(unn_F, 2, max)
    factors_with_unn <- factor_max_unn * 2 > apply(res$F, 2, max)
    if (!any(factors_with_unn)) next

    L_mn <- as.data.frame(poisson2multinom(res)$L)
    colnames(L_mn) <- paste0(gene, ".", colnames(L_mn))
    L_keep <- L_mn[, factors_with_unn, drop = FALSE]
    L_keep <- rownames_to_column(L_keep, "rail_id")

    per_gene_L[[k]] <- column_to_rownames(
      full_join(unfilt_input, L_keep, by = "rail_id"), "rail_id")
  }
  combined <- do.call(cbind, per_gene_L[lengths(per_gene_L) > 0L])
  if (length(combined) == 0L) next
  combined <- rownames_to_column(as.data.frame(combined), "rail_id")
  unfilt   <- unfilt %>% full_join(combined, by = "rail_id")
  message(sprintf("chr%s  features so far: %d", chr, ncol(unfilt) - 1L))
}

write_tsv(unfilt, opt$features_out)
message(sprintf("Wrote %s with %d samples x %d features",
                 opt$features_out, nrow(unfilt), ncol(unfilt) - 1L))

# ---- 2. RIN + productive/unproductive QC ----------------------------------

rin_scores <- read_tsv(opt$analyte, show_col_types = FALSE) %>%
  filter(analytes.rna_integrity_number != "'--") %>%
  transmute(analytes.submitter_id,
             RIN = as.numeric(analytes.rna_integrity_number))

meta_with_rin <- tcga_metadata %>%
  left_join(rin_scores,
             by = c("gdc_cases.samples.portions.analytes.submitter_id"
                    = "analytes.submitter_id")) %>%
  mutate(rail_id_chr = as.character(rail_id)) %>%
  distinct(rail_id_chr, .keep_all = TRUE)

pu_files <- list.files(opt$pu_dir,
                        pattern = "_productive_unproductive\\.tsv$",
                        full.names = TRUE)
combined_pu <- bind_rows(lapply(pu_files, read_tsv, show_col_types = FALSE)) %>%
  group_by(rail_id) %>%
  summarise(productive_counts   = sum(productive_counts),
             unproductive_counts = sum(unproductive_counts),
             .groups = "drop") %>%
  mutate(rail_id_chr      = as.character(rail_id),
          unprod_prod_ratio = unproductive_counts / pmax(productive_counts, 1))

meta_qc <- meta_with_rin %>%
  mutate(pctC            = as.numeric(`%C`),
          avgQ            = as.numeric(avgQ),
          junction_count  = as.numeric(junction_count),
          uniq_mapped_pct = as.numeric(`star.uniquely_mapped_reads_%_both`)) %>%
  left_join(combined_pu %>% dplyr::select(rail_id_chr, unprod_prod_ratio),
             by = "rail_id_chr")

z <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
meta_qc <- meta_qc %>%
  mutate(unprod_z      = z(unprod_prod_ratio),
          pctC_z        = z(pctC),
          RIN_z         = z(RIN),
          avgQ_z        = z(avgQ),
          uniq_mapped_z = z(uniq_mapped_pct))

fail <- meta_qc %>%
  mutate(fail_unprod   = !is.na(unprod_z)      & unprod_z > 3,
          fail_pctC     = !is.na(pctC_z)        & pctC_z   < -3,
          fail_junc     = !is.na(junction_count) & junction_count < 125000,
          fail_rin      = !is.na(RIN)            & RIN < 5,
          fail_avgQ     = !is.na(avgQ_z)         & avgQ_z   < -3,
          fail_uniq_map = !is.na(uniq_mapped_z)  & uniq_mapped_z < -3,
          any_fail = fail_unprod | fail_pctC | fail_junc |
                     fail_rin | fail_avgQ | fail_uniq_map)

removed_rail_ids <- fail$rail_id[fail$any_fail]
message(sprintf("QC: %d samples failing", length(removed_rail_ids)))

# ---- 3. softImpute + flashier ---------------------------------------------

just_factors <- unfilt %>%
  filter(!rail_id %in% c(normal_rail_ids, removed_rail_ids))
to_pca <- column_to_rownames(just_factors, "rail_id")
message(sprintf("Samples after QC + normal removal: %d", nrow(to_pca)))

X <- as.matrix(to_pca)
fit_si    <- softImpute(X, rank.max = min(dim(X)) - 1L, lambda = 0)
X_imputed <- complete(X, fit_si)
temp      <- scale(X_imputed)

out <- flash(temp,
              greedy_Kmax = opt$greedy_kmax,
              ebnm_fn     = list(ebnm_point_exponential, ebnm_point_laplace),
              var_type    = 2,
              backfit     = FALSE)
message(sprintf("flashier: %d factors", out$n_factors))
saveRDS(out, opt$flash_out)
message(sprintf("Wrote %s", opt$flash_out))
