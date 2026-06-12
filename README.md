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
│   └── quetzal-r.yml             # R 4.4 + Bioc 3.20 conda env
├── scripts/
│   ├── gene_level/
│   │   ├── Snakefile             # per-gene Poisson NMF Snakemake workflow
│   │   ├── gene_plots_and_objs.R # fits Poisson NMF + DE per gene
│   │   └── envs/rscript.yml      # symlink-equivalent to ../../environment/quetzal-r.yml
│   └── genome_wide/
│       ├── lf_Snakefile               # leafcutter -> per-gene fastTopics (run once per chr)
│       ├── setting_up_snakemake.sh    # example SLURM dispatcher (edit placeholders before use)
│       ├── tcga_LF_saving.R
│       ├── fasttopics_to_flashier.R   # full FastTopics -> softImpute -> flashier
│       └── envs/fasttopics.yml
├── data/
│   ├── hg38_granges.RDS                  (~24 MB)
│   ├── analyte.tsv                       (~18 MB)
│   ├── all_genes/                        (NOT in repo -- per-chr snaptron junction tables)
│   └── productive_unproductive/          (~8 MB; chr<N>_*.tsv)
└── output/                               (all pipeline outputs land here)
    ├── genome_wide/<chr>/FastTopics_output/<gene>/res.RDS
    └── gene_level/<chr>/<gene>/{res,de_res}.RDS, *.tsv, *.html
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

# create the R environment (R 4.4 + Bioc 3.20 + fastTopics + flashier)
conda env create -f environment/quetzal-r.yml
conda activate quetzal-r
```

`environment/quetzal-r.yml` pins r-base to 4.4 and lets conda solve a
self-consistent Bioconductor 3.20 + CRAN snapshot. The two per-Snakefile
yml files under `scripts/*/envs/` are mirrors of the central yml --
`snakemake --use-conda` builds the env automatically, so you don't need
to run `conda env create` by hand unless you want the env outside of
Snakemake (e.g. for running `fasttopics_to_flashier.R` interactively).

## Methods

Quetzal ships two **independent** methods. They both consume the same
input (per-chromosome snaptron junction tables) and they both fit
Poisson NMFs per gene, but they answer different questions and neither
depends on the other's output. Run one, the other, or both, depending
on what you need.

| Method | Question it answers | Driver | Output location |
| --- | --- | --- | --- |
| **Genome-wide** | "What are the cohort-driving splicing factor programs across all of TCGA?" -- yields one sample x factor matrix (the GSPs) covering every gene at once | `scripts/genome_wide/lf_Snakefile` + `scripts/genome_wide/fasttopics_to_flashier.R` | `output/genome_wide/` |
| **Gene-level** | "For one specific gene, what are its factorised splicing programs?" -- yields per-gene plots and DE results for visual inspection | `scripts/gene_level/Snakefile` (`gene_plots_and_objs.R`) | `output/gene_level/` |

The two pipelines fit independent Poisson NMFs. The genome-wide one
caps at 32 factors per gene and feeds those into a downstream
softImpute + flashier collapse to a single 300-factor matrix. The
gene-level one caps at 10 factors per gene for cleaner per-gene
visualisation and runs its own DE tests on top.

### Genome-wide method

**Step A (Snakemake).** `scripts/genome_wide/lf_Snakefile` consumes one
chromosome's snaptron junction tables and writes a per-gene fastTopics
factorisation. The Snakefile takes the current working directory's
basename as the chromosome label, so the standard invocation is one
Snakemake job per chromosome from a `chr<N>/` subdir of
`scripts/genome_wide/`:

```bash
cd scripts/genome_wide
mkdir chr5 && cp lf_Snakefile chr5/
cd chr5
snakemake -s lf_Snakefile --use-conda --jobs <N>   # add your scheduler flags
```

Outputs land at
`output/genome_wide/<chr>/FastTopics_output/<gene>/res.RDS`.

v0.1 ships no required dispatch wrapper -- parallelise across
chromosomes however your cluster prefers (SLURM array, snakemake
profile, a `for` loop, ...). `setting_up_snakemake.sh` next to
`lf_Snakefile` is one concrete example we used for the manuscript;
edit the two SLURM placeholders at the top of the script (or rewrite
the snakemake invocation for your scheduler) before running.
v1.0 will provide a built-in dispatcher.

**Step B (R script).** Once Step A has populated
`output/genome_wide/<chr>/FastTopics_output/`, collapse every per-gene
result into one sample x feature matrix (keeping only factors that
load mostly on unannotated junctions), QC-filter samples on
RIN / %C / junction count / avgQ / unique-mapped % / unproductive-
productive ratio, fill missing values with `softImpute`, and fit a
300-factor `flashier::flash`:

```bash
conda activate quetzal-r        # only the R env, no snakemake needed
Rscript scripts/genome_wide/fasttopics_to_flashier.R \
    --gene_dir       output/genome_wide \
    --snaptron_root  data/all_genes \
    --tcga_meta      data/tcga_v2_samples.tsv \
    --analyte        data/analyte.tsv \
    --pu_dir         data/productive_unproductive \
    --features_out   output/genome_wide/softimpute_features.tsv \
    --flash_out      output/genome_wide/softimpute_flash_300_qc_filtered.RDS
```

`softimpute_flash_300_qc_filtered.RDS` is the genome-wide factorisation
used by every downstream Quetzal analysis.

External requirement: `data/all_genes/<chr>/snaptron_output/<gene>_snaptron.tsv`
for every gene you want factorised. v0.1 expects this path to exist;
v1.0 will ship a generator/loader for it.

### Gene-level method

`scripts/gene_level/Snakefile` fits its own per-gene Poisson NMF
(low-k cap, for visualisation) and runs differential expression on the
factors:

```bash
cd scripts/gene_level
snakemake --use-conda --jobs <N>   # add your scheduler flags
```

Outputs land at `output/gene_level/<chr>/<gene>/`:
`res.RDS`, `de_res.RDS`, `whole_factor.html`.

This method is independent of the genome-wide method -- it does **not**
consume `output/genome_wide/`. Run it whenever you want a per-gene view
without first running the full genome-wide pipeline.

## License

MIT - see [LICENSE](LICENSE).
