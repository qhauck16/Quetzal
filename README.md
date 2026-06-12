# Quetzal

Per-gene Poisson NMF + genome-wide flashier factorisation of TCGA splice
junctions. Quetzal recovers known splicing-factor mutation programs
(e.g. SF3B1) directly from leafcutter junction counts and surfaces
cohort-specific factors as a single sample x factor matrix that plugs
into downstream alternative-splicing analyses.

This repository is the **v0.1** release used for the manuscript. The
pipelines and reference data here are TCGA-specific; v1.0 will generalise
the sample-metadata assumptions.

## Repository layout

```
Quetzal/
├── environment/
│   └── quetzal-r.yml             # pinned R env (R 4.2.3)
├── scripts/
│   ├── gene_level/
│   │   ├── Snakefile             # per-gene Poisson NMF Snakemake workflow
│   │   ├── gene_plots_and_objs.R # fits Poisson NMF + DE per gene
│   │   ├── cancer_specific_factors.R
│   │   └── envs/rscript.yml      # symlink-equivalent to ../../environment/quetzal-r.yml
│   └── genome_wide/
│       ├── lf_Snakefile               # leafcutter -> per-gene fastTopics (run once per chr)
│       ├── setting_up_snakemake.sh    # example SLURM dispatcher (edit placeholders before use)
│       ├── tcga_LF_saving.R
│       ├── fasttopics_to_flashier.R   # full FastTopics -> softImpute -> flashier
│       └── envs/fasttopics.yml
└── data/
    ├── hg38_granges.RDS                  (~24 MB)
    ├── analyte.tsv                       (~18 MB)
    └── productive_unproductive/          (~8 MB; chr<N>_*.tsv)
```

### Data not shipped with the repo

The Snaptron TCGA-v2 sample metadata file is **not** included in this
repository (~65 MB, distributed by Snaptron upstream). Download it
once into `data/` before running any of the pipelines:

```bash
curl -L -o data/tcga_v2_samples.tsv \
     https://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv
```

## Install

```bash
git clone <repo-url>
cd Quetzal

# create the R environment (R 4.2.3 + flashier 1.0.58 + fastTopics 0.6-192 + ...)
conda env create -f environment/quetzal-r.yml
conda activate quetzal-r
```

`environment/quetzal-r.yml` pins every R package to the exact version
used at the v0.1 freeze. If you only want a specific stage you can skip
either of the snakefile envs -- they alias the central yml.

## End-to-end run (TCGA)

The full pipeline has two stages with independent Snakemake workflows.

### 1. Per-chromosome junction files (v0.1 assumes Snaptron-format) -> fastTopics (genome-wide)

`scripts/genome_wide/lf_Snakefile` consumes one chromosome's snaptron
junction tables (produced upstream from recount3 / Snaptron) and writes
a fastTopics factorisation per gene under
`scripts/genome_wide/<chr>/FastTopics_output/<gene>/res.RDS`.

The Snakefile takes the current working directory's basename as the
chromosome label, so the standard invocation is one Snakemake job per
chromosome from a `chr<N>/` subdir of `scripts/genome_wide/`:

```bash
cd scripts/genome_wide
mkdir chr5 && cp lf_Snakefile chr5/
cd chr5
snakemake -s lf_Snakefile --use-conda --jobs <N>   # add your scheduler flags
```

External requirement: `data/all_genes/<chr>/snaptron_output/<gene>_snaptron.tsv`
for every gene you want factorised. v0.1 expects this path to exist;
v1.0 will ship a generator/loader for it.

v0.1 ships no required dispatch wrapper -- parallelise across
chromosomes however your cluster prefers (SLURM array, snakemake
profile, a `for` loop, ...). `setting_up_snakemake.sh` next to
`lf_Snakefile` is one concrete example we used for the manuscript;
edit the two SLURM placeholders at the top of the script (or rewrite
the snakemake invocation for your scheduler) before running.
v1.0 will provide a built-in dispatcher.

### 2. Per-gene Poisson NMF + DE (gene-level)

`scripts/gene_level/Snakefile` reads each gene's fastTopics output from
step 1 and runs `gene_plots_and_objs.R` (Poisson NMF + DE) then
`cancer_specific_factors.R` (beta-binomial cohort test):

```bash
cd scripts/gene_level
snakemake --use-conda --jobs <N>   # add your scheduler flags
```

### 3. Gene-level -> genome-wide flashier

`scripts/genome_wide/fasttopics_to_flashier.R` collapses every per-gene
result into one sample x feature matrix (only factors that load mostly
on unannotated junctions), QC-filters samples on RIN / %C / junction
count / avgQ / unique-mapped % / unproductive-productive ratio, fills
missing values with `softImpute`, and fits a 300-factor `flashier::flash`:

```bash
Rscript scripts/genome_wide/fasttopics_to_flashier.R \
    --gene_dir   ../260116_filters \
    --snaptron_root ../../all_genes \
    --tcga_meta  data/tcga_v2_samples.tsv \
    --analyte    data/analyte.tsv \
    --pu_dir     data/productive_unproductive \
    --features_out softimpute_features.tsv \
    --flash_out  softimpute_flash_300_qc_filtered.RDS
```

The final `softimpute_flash_300_qc_filtered.RDS` is the genome-wide
factorisation used by every downstream Quetzal analysis.

## License

MIT - see [LICENSE](LICENSE).
