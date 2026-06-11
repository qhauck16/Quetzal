.libPaths('/scratch/midway2/qhauck/conda_env/rstudio-server/lib/R/library/')

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(bit64)
  library(aod)            # betabin (beta-binomial GLM)
})

# ---- Args --------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
new_res_path        <- args[1]
de_res_path         <- args[2]
tcga_metadata_path    <- args[3]
gtex_metadata_path    <- args[4]
covariate_table_path  <- args[5]
psi_dir               <- args[6]
lfc_thresh            <- as.numeric(args[7])
abs_thresh            <- as.numeric(args[8])
sval_thresh           <- as.numeric(args[9])
factor_corr_thresh    <- as.numeric(args[10])
rare_factor_thresh    <- as.numeric(args[11])
psi_pval_thresh       <- as.numeric(args[12])
min_num_in_cancer     <- as.integer(args[13])
chr                   <- args[14]
gene_name             <- args[15]
gtex_dat_dir          <- args[16]
# Optional output_base prefix (args[17]); when given, OUTPUT writes go
# under it. `chr` itself stays raw because lines 185/301/302 use it as
# part of INPUT paths (gtex_dat_dir/<chr>, psi_dir/<source>/<chr>).
output_base           <- if (length(args) >= 17 && nzchar(args[17])) args[17] else "."

# ---- Output paths ------------------------------------------------
out_dir       <- file.path(output_base, chr, gene_name)
results_path  <- file.path(out_dir, "results.tsv")
junction_path <- file.path(out_dir, "junction_results.tsv")
n_bb_path     <- file.path(out_dir, "n_bb_tests.tsv")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Always create the declared outputs even when there are no hits / no
# GTEx match / etc. n_bb_tests is written here with a 0 default; we
# overwrite at the end of the BB loop with the actual count so
# gene_cancer_type_hits.Rmd can apply the right Bonferroni denominator.
file.create(results_path)
file.create(junction_path)
n_bb_tests <- 0L
write.table(data.frame(n_bb_tests = n_bb_tests), n_bb_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

# Empty new_res input -> nothing to do
if (file.size(new_res_path) == 0) quit(save = "no", status = 0)

# ---- Load --------------------------------------------------------
new_res        <- readRDS(new_res_path)
de_res         <- readRDS(de_res_path)
tcga_metadata  <- fread(tcga_metadata_path, data.table = FALSE)
gtex_metadata  <- fread(gtex_metadata_path, data.table = FALSE)

dat <- fread(file.path(out_dir, "dat.tsv"), data.table = FALSE) %>%
  column_to_rownames("V1") %>% as.matrix()

# ---- TCGA -> GTEx tissue lookup ----------------------------------
tcga_to_gtex <- c(
  "TCGA-ACC"   = "ADRENAL_GLAND", "TCGA-BLCA"  = "BLADDER",
  "TCGA-BRCA"  = "BREAST",        "TCGA-CESC"  = "CERVIX_UTERI",
  "TCGA-COAD"  = "COLON",
  "TCGA-ESCA"  = "ESOPHAGUS",     "TCGA-GBM"   = "BRAIN", "TCGA-KICH"  = "KIDNEY",
  "TCGA-KIRC"  = "KIDNEY",        "TCGA-KIRP"  = "KIDNEY",
  "TCGA-LAML"  = "BLOOD",         "TCGA-LGG"   = "BRAIN",
  "TCGA-LIHC"  = "LIVER",         "TCGA-LUAD"  = "LUNG",
  "TCGA-LUSC"  = "LUNG",          "TCGA-MESO"  = "LUNG",
  "TCGA-OV"    = "OVARY",         "TCGA-PAAD"  = "PANCREAS", "TCGA-PRAD"  = "PROSTATE",
  "TCGA-READ"  = "COLON",         "TCGA-SARC"  = "MUSCLE",
  "TCGA-SKCM"  = "SKIN",          "TCGA-STAD"  = "STOMACH",
  "TCGA-TGCT"  = "TESTIS",        "TCGA-THCA"  = "THYROID",
  "TCGA-UCEC"  = "UTERUS",        "TCGA-UCS"   = "UTERUS"
)

# -----------------------------------------------------------------
# Covariate handling.
#
# Covariates are read from the standardized table built by
# tcga/covariate_table.R. RIN is dropped here because it is missing for
# whole cohorts (TCGA-LAML, etc.) and is perfectly collinear with cohort
# identity for any imputation scheme; we'd rather lose the RIN signal
# than lose those cohorts to a downstream complete.cases() filter. We do
# no imputation -- samples with any remaining NA covariate are dropped
# below.
# -----------------------------------------------------------------
pre_model <- readRDS(covariate_table_path) %>%
  dplyr::select(-analytes.rna_integrity_number) %>%
  mutate(rail_id = as.numeric(rail_id))

# -----------------------------------------------------------------
# Match NMF L to covariates, drop incomplete rows + tiny cancer types
# -----------------------------------------------------------------
nmf_results <- as.data.frame(new_res$L) %>%
  rownames_to_column("rail_id") %>%
  mutate(rail_id = as.numeric(rail_id))

full_covar <- pre_model %>%
  filter(rail_id %in% nmf_results$rail_id) %>%
  filter(complete.cases(.)) %>%
  group_by(gdc_cases.project.project_id) %>%
  filter(n() >= min_num_in_cancer) %>%
  ungroup()

# Drop columns that are constant across what's left
full_covar <- full_covar[, sapply(full_covar,
                                  function(x) length(unique(x)) > 1)]

rail_ids   <- full_covar$rail_id
full_covar <- dplyr::select(full_covar, -rail_id)

# -----------------------------------------------------------------
# Build covariate model matrix; remove colinear columns via QR
# -----------------------------------------------------------------
covar_matrix <- model.matrix(~ ., data = full_covar)
qr_d  <- qr(covar_matrix)
nc    <- ncol(covar_matrix)
if (qr_d$rank < nc) {
  drop_idx  <- qr_d$pivot[(qr_d$rank + 1):nc]
  drop_name <- colnames(covar_matrix)[drop_idx]
  message("Removing colinear columns: ", paste(drop_name, collapse = ", "))
  covar_matrix_final <- covar_matrix[, -drop_idx]
} else {
  message("No colinear columns detected.")
  covar_matrix_final <- covar_matrix
}

# -----------------------------------------------------------------
# Per-sample frame: one row per rail_id with covariates + factors
# -----------------------------------------------------------------
pre_model_filtered <- as.data.frame(covar_matrix_final) %>%
  mutate(rail_id = rail_ids) %>%
  left_join(nmf_results, by = "rail_id") %>%
  column_to_rownames("rail_id")

factor_cols <- grep("^k\\d+$",                names(pre_model_filtered), value = TRUE)
cancer_cols <- grep("project.project_id",     names(pre_model_filtered), value = TRUE)

# Only test "rare" factors (mean L below threshold)
tcga_factor_means <- colMeans(as.data.frame(new_res$L))
rare_factors      <- names(which(tcga_factor_means < rare_factor_thresh))

# -----------------------------------------------------------------
# Linear model: each rare factor x each cancer type
# -----------------------------------------------------------------
hits        <- list()
coef_tables <- list()
hit_idx     <- 0L

for (fact in rare_factors) {
  this_k <- pre_model_filtered[[fact]]
  for (cancer_col in cancer_cols) {
    this_cancer <- pre_model_filtered[[cancer_col]]

    just_this_cancer <- pre_model_filtered %>%
      dplyr::select(-all_of(cancer_cols), -all_of(factor_cols)) %>%
      mutate(cancer_statusTRUE = this_cancer,
             k                 = this_k)

    fit <- lm(k ~ ., data = just_this_cancer)
    co  <- summary(fit)$coefficients
    if (!"cancer_statusTRUE" %in% rownames(co)) next

    pval <- co["cancer_statusTRUE", 4]
    beta <- co["cancer_statusTRUE", 1]

    if (pval < factor_corr_thresh && beta > 0) {
      hit_idx <- hit_idx + 1L
      hits[[hit_idx]] <- list(
        factor      = fact,
        cancer_type = str_split_i(cancer_col, "_id", 2),
        cancer_pval = pval,
        cancer_beta = beta
      )
      coef_tables[[hit_idx]] <- co
    }
  }
}

if (length(hits) == 0) quit(save = "no", status = 0)

# -----------------------------------------------------------------
# GTEx data setup
# -----------------------------------------------------------------
gtex_dat_path <- file.path(gtex_dat_dir, chr,
                           "FastTopics_output", gene_name, "dat.tsv")
if (!file.exists(gtex_dat_path) || file.size(gtex_dat_path) < 100) {
  quit(save = "no", status = 0)
}

gtex_dat <- fread(gtex_dat_path, data.table = FALSE) %>%
  column_to_rownames("rownames")
gtex_dat <- gtex_dat[, colnames(gtex_dat) %in% colnames(dat)]
for (col in setdiff(colnames(dat), colnames(gtex_dat))) gtex_dat[[col]] <- 0
gtex_dat <- as.matrix(gtex_dat[, c(colnames(dat),
                                   setdiff(colnames(gtex_dat),
                                           colnames(dat)))])

# -----------------------------------------------------------------
# Beta-binomial GLM helper
#
# Returns NULL if the model can't be fit. Always produces per-group
# PSI estimates (logit-1 of intercept and intercept + studyTCGA).
# -----------------------------------------------------------------
fit_bb <- function(cluster_test) {
  if (length(unique(cluster_test$study)) < 2) return(NULL)
  cluster_test$study <- factor(cluster_test$study,
                               levels = c("GTEx", "TCGA"))

  res <- tryCatch(
    suppressWarnings(aod::betabin(
      cbind(junc_count, cluster_sum - junc_count) ~ study,
      random = ~ 1,
      data   = cluster_test
    )),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  co <- summary(res)@Coef
  if (!"studyTCGA" %in% rownames(co)) return(NULL)

  # Locate the p-value column robustly across aod versions
  p_col <- grep("^Pr|p.value|pvalue", colnames(co),
                ignore.case = TRUE, value = TRUE)[1]
  if (is.na(p_col)) return(NULL)

  est_int  <- co["(Intercept)", "Estimate"]
  est_tcga <- co["studyTCGA",   "Estimate"]
  list(
    beta     = est_tcga,
    se       = co["studyTCGA", "Std. Error"],
    pvalue   = co["studyTCGA", p_col],
    psi_gtex = plogis(est_int),
    psi_tcga = plogis(est_int + est_tcga)
  )
}

# -----------------------------------------------------------------
# Per-junction boxplot helper (unchanged behaviour, factored out)
# -----------------------------------------------------------------
boxplot_one_junc <- function(junc, gtex_norm_dat, tcga_norm_dat,
                             gtex_label, cancer_label, junc_pval) {
  gtex_idx <- which(grepl(junc, colnames(gtex_norm_dat)))
  tcga_idx <- which(grepl(junc, colnames(tcga_norm_dat)))
  out <- list()
  for (k in seq_along(gtex_idx)) {
    to_plot <- data.frame(
      junc      = c(gtex_norm_dat[, gtex_idx[k]],
                    tcga_norm_dat[, tcga_idx[k]]),
      cohort_id = c(rep(gtex_label,   nrow(gtex_norm_dat)),
                    rep(cancer_label, nrow(tcga_norm_dat)))
    )
    p_lab <- if (!is.na(junc_pval))
      sprintf("p = %.3g", junc_pval) else "p = NA"
    out[[length(out) + 1]] <- ggplot(to_plot, aes(x = cohort_id, y = junc)) +
      geom_boxplot() +
      labs(title = paste0(gene_name, ":", junc),
           y     = "Leafcutter PSI") +
      theme_minimal() +
      annotate("text",
               x     = 1.5,
               y     = max(to_plot$junc, na.rm = TRUE) * 1.05,
               label = p_lab, size = 4, hjust = 0.5)
  }
  out
}

# -----------------------------------------------------------------
# Per-hit junction analysis + outputs
# -----------------------------------------------------------------
all_junction_rows <- list()
results_rows      <- list()

for (i in seq_along(hits)) {
  h           <- hits[[i]]
  fact        <- h$factor
  cancer_type <- h$cancer_type
  if (!cancer_type %in% names(tcga_to_gtex)) next
  gtex_study  <- unname(tcga_to_gtex[cancer_type])

  # subset count matrices to this cancer type / GTEx tissue
  tcga_rids <- filter(tcga_metadata,
                      gdc_cases.project.project_id == cancer_type)$rail_id
  gtex_rids <- filter(gtex_metadata, study == gtex_study)$rail_id
  dat_cancer      <- dat[rownames(dat) %in% tcga_rids, , drop = FALSE]
  dat_gtex_tissue <- gtex_dat[rownames(gtex_dat) %in% gtex_rids, , drop = FALSE]

  # candidate junctions for this factor
  fact_idx <- as.integer(sub("^k", "", fact))
  svals    <- de_res$svalue;   svals[is.na(svals)] <- 1
  lfcs     <- de_res$postmean; lfcs [is.na(lfcs)]  <- 0
  junctions_to_compare <- colnames(dat)[
    svals[, fact_idx] < sval_thresh &
    lfcs [, fact_idx] > lfc_thresh &
    new_res$F[, fact_idx] > abs_thresh
  ]
  if (length(junctions_to_compare) == 0) next

  # PSI tables for graphing / reporting
  gtex_psi_path <- file.path(psi_dir, "gtex", chr, gene_name, "psi_table.tsv")
  tcga_psi_path <- file.path(psi_dir, "tcga", chr, gene_name, "psi_table.tsv")
  if (!file.exists(gtex_psi_path) || file.size(gtex_psi_path) < 100) next
  if (!file.exists(tcga_psi_path) || file.size(tcga_psi_path) < 100) next
  gtex_norm_dat <- fread(gtex_psi_path, data.table = FALSE) %>%
    filter(V1 %in% rownames(dat_gtex_tissue)) %>% column_to_rownames("V1")
  tcga_norm_dat <- fread(tcga_psi_path, data.table = FALSE) %>%
    filter(V1 %in% rownames(dat_cancer)) %>% column_to_rownames("V1")

  # ---- BB-test every candidate junction ----
  junc_rows <- list()
  for (junc in junctions_to_compare) {
    tcga_clust <- str_split_i(names(tcga_norm_dat)[grepl(junc, names(tcga_norm_dat))], ":", 1)
    gtex_clust <- str_split_i(names(gtex_norm_dat)[grepl(junc, names(gtex_norm_dat))], ":", 1)
    clust <- intersect(tcga_clust, gtex_clust)[1]
    if (is.na(clust)) next

    tcga_clust_juncs <- str_split_i(names(tcga_norm_dat)[grepl(clust, names(tcga_norm_dat))], ":", 2)
    gtex_clust_juncs <- str_split_i(names(gtex_norm_dat)[grepl(clust, names(gtex_norm_dat))], ":", 2)

    tcga_vals <- dat[rownames(dat) %in% rownames(tcga_norm_dat),
                     junc, drop = FALSE]
    tcga_sums <- rowSums(dat[rownames(dat) %in% rownames(tcga_norm_dat),
                             colnames(dat) %in% tcga_clust_juncs,
                             drop = FALSE])
    gtex_vals <- gtex_dat[rownames(gtex_dat) %in% rownames(gtex_norm_dat),
                          junc, drop = FALSE]
    gtex_sums <- rowSums(gtex_dat[rownames(gtex_dat) %in% rownames(gtex_norm_dat),
                                  colnames(gtex_dat) %in% gtex_clust_juncs,
                                  drop = FALSE])

    cluster_test <- data.frame(
      junc_count  = c(tcga_vals, gtex_vals),
      cluster_sum = c(tcga_sums, gtex_sums),
      study       = c(rep("TCGA", length(tcga_vals)),
                      rep("GTEx", length(gtex_vals))),
      rail_id     = c(rownames(tcga_vals), rownames(gtex_vals)),
      stringsAsFactors = FALSE
    ) %>% filter(cluster_sum > 0)

    bb <- fit_bb(cluster_test)
    if (is.null(bb)) next
    n_bb_tests <- n_bb_tests + 1L     # count every successful BB fit

    psi_id <- paste0(clust, ":", junc)
    median_psi_tcga <- median(tcga_norm_dat[, psi_id], na.rm = TRUE)
    median_psi_gtex <- median(gtex_norm_dat[, psi_id], na.rm = TRUE)

    junc_rows[[length(junc_rows) + 1]] <- data.frame(
      gene_name             = gene_name,
      chr                   = chr,
      cancer_type           = cancer_type,
      gtex_study            = gtex_study,
      factor                = fact,
      factor_cancer_pvalue  = h$cancer_pval,
      factor_cancer_beta    = h$cancer_beta,
      junction_id           = junc,
      cluster_id            = clust,
      n_tcga                = sum(cluster_test$study == "TCGA"),
      n_gtex                = sum(cluster_test$study == "GTEx"),
      bb_beta               = bb$beta,
      bb_se                 = bb$se,
      bb_pvalue             = bb$pvalue,
      bb_psi_gtex           = bb$psi_gtex,
      bb_psi_tcga           = bb$psi_tcga,
      median_psi_tcga       = median_psi_tcga,
      median_psi_gtex       = median_psi_gtex,
      median_psi_dif        = median_psi_tcga - median_psi_gtex,
      median_psi_fc         = (median_psi_tcga + 1e-4) /
                              (median_psi_gtex + 1e-4),
      stringsAsFactors      = FALSE
    )
  }

  if (length(junc_rows) == 0) next
  junc_df <- bind_rows(junc_rows)

  # require at least one positive-direction significant junction at the
  # uncorrected per-junction threshold (matches old behaviour)
  if (sum(junc_df$bb_pvalue < psi_pval_thresh & junc_df$bb_beta > 0,
          na.rm = TRUE) == 0) next

  # ---- per-hit boxplot PDF (boxplots overlay BB p-values) ----
  pval_lookup <- setNames(junc_df$bb_pvalue, junc_df$junction_id)
  plots <- list()
  for (junc in junctions_to_compare) {
    plots <- c(plots, boxplot_one_junc(junc,
                                       gtex_norm_dat, tcga_norm_dat,
                                       gtex_study, cancer_type,
                                       pval_lookup[junc]))
  }
  if (length(plots) > 0) {
    combined_plot <- patchwork::wrap_plots(plots, ncol = 3) +
      patchwork::plot_annotation(
        title = paste0("Factor ", fact, " is correlated with ", cancer_type,
                       " with pval ", signif(h$cancer_pval, 2)),
        theme = theme(plot.title = element_text(hjust = 0.5,
                                                size = 16, face = "bold"))
      )
    ggplot2::ggsave(file.path(out_dir,
                              paste0(fact, "_", cancer_type, ".pdf")),
                    combined_plot, width = 12, height = 8)
  }

  # ---- per-hit coef table TSV (linear model summary) ----
  write.table(coef_tables[[i]],
              file.path(out_dir, paste0(fact, "_", cancer_type, ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)

  # ---- per-hit junction-level TSV ----
  fwrite(junc_df,
         file.path(out_dir,
                   paste0(fact, "_", cancer_type, "_junctions.tsv")),
         sep = "\t")

  # ---- one row in the gene-level results.tsv (legacy schema) ----
  results_rows[[length(results_rows) + 1]] <- data.frame(
    gene_name              = gene_name,
    cancer_type            = cancer_type,
    factor                 = fact,
    factor_coefficient     = h$cancer_beta,
    pvalue                 = signif(h$cancer_pval, 2),
    num_sig_juncs_in_factor = sum(junc_df$bb_pvalue < psi_pval_thresh,
                                  na.rm = TRUE),
    max_delta_psi          = max(junc_df$median_psi_dif, na.rm = TRUE),
    max_psi_fc             = max(junc_df$median_psi_fc,  na.rm = TRUE),
    chr                    = chr,
    stringsAsFactors       = FALSE
  )

  all_junction_rows[[length(all_junction_rows) + 1]] <- junc_df
}

# Append legacy-schema rows to results.tsv
if (length(results_rows) > 0) {
  fwrite(bind_rows(results_rows), results_path, sep = "\t",
         col.names = FALSE, append = TRUE)
}

# Write per-gene junction-level rollup (one file across all hits)
if (length(all_junction_rows) > 0) {
  fwrite(bind_rows(all_junction_rows), junction_path, sep = "\t")
}

# Per-gene BB test count for downstream Bonferroni. Counts every BB fit
# that succeeded (including those in (factor, cohort) hits that ended up
# being dropped for having zero TCGA-up significant junctions).
write.table(data.frame(n_bb_tests = n_bb_tests), n_bb_path,
            sep = "\t", quote = FALSE, row.names = FALSE)
