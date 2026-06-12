# Quetzal

Per-gene Poisson NMF and genome-wide flashier factorisation of RNA-seq
splice junctions. Quetzal turns junction-count matrices into
sample × factor matrices that surface cohort-driving splicing programs
(e.g. SF3B1 mutation signatures) on their own.

v1.0 generalises Quetzal beyond TCGA+Snaptron: it runs on user-supplied
per-gene junction matrices in addition to snaptron's TSVs, drops every
TCGA-specific QC step, and lets you opt in (or out) of normal-sample
filtering and structure-plot grouping via a single configfile. v0.1
(the manuscript build) is frozen at the `v0.1.0` tag —
`git checkout v0.1.0` if you need that exact pipeline.

## Quickstart

```bash
git clone https://github.com/qhauck16/Quetzal.git
cd Quetzal

# (snaptron input only) fetch the TCGA-v2 sample manifest -- not shipped
curl -L -o data/tcga_v2_samples.tsv \
     https://snaptron.cs.jhu.edu/data/tcgav2/samples.tsv

# (snaptron input only) drop per-chr snaptron tables under data/all_genes/<chr>/
# (user-supplied gene_matrix input goes under the same dir; see below)

# pick your mode + run. Snakemake builds the conda env via --use-conda.
cd scripts/genome_wide      # or scripts/gene_level
snakemake --configfile ../../config/default_config.yaml --use-conda --jobs <N>
```

## Two modes

Run one, the other, or both. Both modes share canonical per-gene
intermediates so the input-format and filter choices look the same in
either mode.

| Mode | Question it answers | Driver | Final output |
| --- | --- | --- | --- |
| **gene_level** | "For one specific gene, what are its factorised splicing programs?" -- per-gene plots and DE results for visual inspection | `scripts/gene_level/Snakefile` | `output/gene_level/<chr>/<gene>/whole_factor.html` (+ `res.RDS`, `de_res.RDS`, `de_factor.html`) |
| **genome_wide** | "What are the cohort-driving splicing programs across all genes?" -- one sample × factor matrix that covers the entire input | `scripts/genome_wide/Snakefile` | `output/genome_wide/softimpute_flash.RDS` (+ per-gene `res.RDS`, `unann_factors.tsv`) |

Both fit independent Poisson NMFs per gene. The gene_level mode caps
`k` at 10 by default for cleaner visualisation; the genome_wide mode
caps higher (32 is typical) and feeds those factors into a downstream
softImpute + flashier collapse. Toggle via the configfile.

## Two input formats

Picked at run time via the `input_format` key in your config:

```yaml
input_format: snaptron       # OR: gene_matrix
input_dir:    data/all_genes # both formats expect per-chromosome subdirs
```

### `snaptron` (default)

Per-gene snaptron TSV at
`<input_dir>/<chr>/snaptron_output/<gene>_snaptron.tsv`. The required
columns are the standard snaptron set:

```
chromosome  start  end  strand  samples  samples_count  annotated  ...
```

The `samples` column is the packed `,railid1:count1,railid2:count2,...`
string. `annotated` (0/1) is taken straight from snaptron and is the
basis for `extract_unann_factors.R`'s annotated/unannotated split.

### `gene_matrix` (user-supplied)

Per-gene wide TSV at `<input_dir>/<chr>/<gene>.tsv`:

```
junction_id              S1   S2   S3   ...
chr5:106932115-106934550:+  4   13    0   ...
chr5:106932115-106934700:+  0    2    7   ...
```

* Column 1: `junction_id` formatted exactly as `chr:start-end:strand`.
* Remaining columns: sample IDs you choose (numeric counts in each cell).

`classify_junctions.R` runs automatically after ingest and fills in the
`annotated` (0/1) flag for each junction by comparing its coords to the
gencode model. That rule is skipped entirely when `filter_unannotated`
is `false` — you can run gene_matrix input through genome_wide mode
without ever touching gencode-derived annotation.

## Sample metadata (optional)

Quetzal needs **nothing** from a sample metadata file in the default
config. Provide it only if you want one of the two opt-in features:

| Feature | Mode | Required config | Required metadata columns |
| --- | --- | --- | --- |
| Drop normal/control samples upfront | both | `exclude_normals: true` + `normal_filter.column` + `normal_filter.pattern` | `sample_id_column`, `normal_filter.column` |
| Group samples in the structure plot | gene_level | `structure_plot_grouping_column` | `sample_id_column`, `structure_plot_grouping_column` |

`sample_id_column` (default `rail_id`) must match either the snaptron
rail_ids inside the packed `samples` column (snaptron input) or the
column headers of your wide TSV (gene_matrix input).

### Shipped example: `data/example_tcga_metadata.tsv`

A 348 KB slim version of the Snaptron TCGA-v2 manifest with just the
three columns Quetzal actually consumes for TCGA cohorts:

```
rail_id  gdc_cases.project.project_id  cgc_sample_sample_type
106797   TCGA-ACC                      Primary Tumor
110230   TCGA-ACC                      Primary Tumor
...
```

The matching config snippet to reproduce v0.1 behaviour on TCGA:

```yaml
sample_metadata: data/example_tcga_metadata.tsv
sample_id_column: rail_id
exclude_normals: true
normal_filter:
  column: cgc_sample_sample_type
  pattern: Normal
# gene_level only: group structure plot by TCGA cancer cohort
structure_plot_grouping_column: gdc_cases.project.project_id
```

Drop in your own manifest in the same shape for non-TCGA cohorts —
swap `cgc_sample_sample_type` / `gdc_cases.project.project_id` for
whatever your dataset calls those concepts.

## Configfile

A single `config/default_config.yaml` is the source of truth for every
knob. Copy it, edit, point `snakemake --configfile your_config.yaml` at
it. The shipped defaults:

```yaml
# input
input_format: snaptron        # or gene_matrix
input_dir:    data/all_genes
output_dir:   output
gencode:      data/gencode_v46_granges.RDS

# per-gene NMF (both modes)
max_factors:                       10      # gene_level default; bump for genome_wide (~32)
variance_explained:                0.99    # PCA cumulative-variance elbow
min_samples_per_junc:              10
min_reads_per_sample_per_cluster:  5
min_clust_read_count_avg:          0
gene_range_bound:                  2000

# optional sample-level filtering
sample_metadata:    null
sample_id_column:   rail_id
exclude_normals:    false
normal_filter:
  column:  null
  pattern: null

# gene_level structure-plot grouping (optional)
structure_plot_grouping_column: null

# genome_wide
filter_unannotated:           true   # keep only factors loading on unannotated junctions
unann_factor_loading_ratio:   2.0
sample_fraction:              0.8    # min coverage to include a gene's factors
greedy_kmax:                  300    # flashier::flash greedy_Kmax
```

## Running the pipelines

Snakemake builds the conda env from `scripts/<mode>/envs/rscript.yml`
(a mirror of `environment/quetzal-r.yml`) on first invocation. Add
`--executor slurm --default-resources slurm_partition=... slurm_account=...`
for cluster fan-out.

```bash
# gene_level
cd scripts/gene_level
snakemake --configfile ../../config/default_config.yaml --use-conda --jobs <N>

# genome_wide
cd scripts/genome_wide
snakemake --configfile ../../config/default_config.yaml --use-conda --jobs <N>
```

Each per-gene rule writes a `list(skipped = TRUE, reason = "...")`
stub to its output when filters knock the gene out, so the DAG never
breaks — failed genes just propagate as empty rows in the final
flashier matrix.

### v0.1-equivalent config (for reproducing manuscript outputs)

```yaml
input_format:   snaptron
input_dir:      data/all_genes
output_dir:     output
gencode:        data/gencode_v46_granges.RDS

# v0.1's elbow_cutoff = 0.01 == v1.0's variance_explained = 0.99
variance_explained:                0.99
max_factors:                       10
min_samples_per_junc:              10
min_reads_per_sample_per_cluster:  5
gene_range_bound:                  2000

# v0.1 silently dropped TCGA normals
sample_metadata:  data/example_tcga_metadata.tsv
sample_id_column: rail_id
exclude_normals:  true
normal_filter:
  column:  cgc_sample_sample_type
  pattern: Normal

# v0.1 grouped the structure plot by TCGA cancer cohort
structure_plot_grouping_column: gdc_cases.project.project_id
```

## Repository layout

```
Quetzal/
├── README.md
├── LICENSE
├── config/
│   └── default_config.yaml
├── environment/
│   └── quetzal-r.yml                  # R 4.4 + Bioc 3.20 conda env
├── scripts/
│   ├── common/                        # shared adapters + NMF + per-gene unann extract
│   │   ├── _canonical.R               # post-parse pipeline (shared by both ingests)
│   │   ├── ingest_snaptron.R
│   │   ├── ingest_gene_matrix.R
│   │   ├── classify_junctions.R
│   │   ├── fit_pnmf.R
│   │   └── extract_unann_factors.R
│   ├── gene_level/
│   │   ├── Snakefile
│   │   ├── de_analysis.R
│   │   ├── make_plots.R
│   │   └── envs/rscript.yml           # mirror of environment/quetzal-r.yml
│   └── genome_wide/
│       ├── Snakefile
│       ├── softimpute_flash.R
│       └── envs/rscript.yml
├── data/
│   ├── gencode_v46_granges.RDS               # default gencode (24 MB)
│   ├── example_tcga_metadata.tsv      # 348 KB slim TCGA manifest
│   └── all_genes/                     # per-chr inputs (NOT shipped; user-provided)
└── output/                            # snakemake target dir
    ├── gene_level/<chr>/<gene>/{matrix,res,de_res}.RDS, *.html
    └── genome_wide/<chr>/<gene>/{matrix,res}.RDS, unann_factors.tsv,
        softimpute_flash.RDS
```

## License

MIT — see [LICENSE](LICENSE).
