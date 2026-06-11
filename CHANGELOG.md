# Changelog

## v0.1.0 - 2026-06-11

Initial release accompanying the manuscript-stage run on TCGA.

### What's included
- `scripts/genome_wide/` — leafcutter junction → per-gene fastTopics
  pipeline (`lf_Snakefile`, `setting_up_snakemake.sh`, `tcga_LF_saving.R`)
  and the `fasttopics_to_flashier.R` end-to-end script that turns the
  per-gene FastTopics outputs into the final 300-factor flashier object
  with QC filtering (RIN, %C, junction count, avgQ, unique-mapped %,
  unproductive/productive-junction ratio).
- `scripts/gene_level/` — `Snakefile` running `gene_plots_and_objs.R`
  and `cancer_specific_factors.R` per gene, plus a beta-binomial test
  for cancer-type-specific factor loadings.
- `environment/quetzal-r.yml` — conda env pinning the exact R / package
  versions used at the v0.1 freeze.
- `data/` — TCGA v2 sample metadata, hg38 GRanges object, RIN-score
  `analyte.tsv`, and the per-chromosome NMD productive/unproductive
  read counts used by the QC filter.

### Known v0.1 limitations
The pipeline currently hard-codes assumptions specific to the TCGA v2
sample list used for the paper:
- Sample metadata column names (`gdc_cases.*`, `cgc_sample_*`,
  `star.uniquely_mapped_reads_*`) are GDC/recount-specific.
- The "normal" sample exclusion looks for the literal string `Normal`
  in `cgc_sample_sample_type`.
- The per-gene snaptron junction tables (under `all_genes/`) are NOT
  shipped here -- they're TCGA-derived ~10 GB of data; v0.1 expects
  them at the path passed via `--snaptron_root`.

### Planned for v1.0
- Make sample metadata column names configurable (or auto-detected).
- Distribute or auto-fetch the per-gene snaptron junction tables.
- Add a small synthetic example dataset that runs end-to-end in <1 min
  so the install can be smoke-tested.
- Replace the per-script ad-hoc CLI parsing with a single config YAML.
- CI: run the smoke test on every PR to main.
