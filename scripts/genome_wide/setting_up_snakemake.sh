#!/usr/bin/env bash
# Quetzal v0.1 -- per-chromosome leafcutter -> fastTopics dispatcher.
#
# Run from <repo>/scripts/genome_wide/. For every chr<N>/ subdir present
# under <repo>/data/all_genes/, materialise <repo>/scripts/genome_wide/<chr>/,
# drop a copy of lf_Snakefile in it, and submit one Snakemake job
# (SLURM, --use-conda) restricted to that chromosome.

set -euo pipefail

ALL_GENES_DIR="../../data/all_genes"

if [ ! -d "$ALL_GENES_DIR" ]; then
    echo "error: $ALL_GENES_DIR does not exist." >&2
    echo "Place per-chr snaptron junction tables at" >&2
    echo "  data/all_genes/<chr>/snaptron_output/<gene>_snaptron.tsv" >&2
    exit 1
fi

for chr_dir in "$ALL_GENES_DIR"/chr*/; do
    chr=$(basename "$chr_dir")
    mkdir -p "$chr"
    cp lf_Snakefile "$chr/"

    cd "$chr"
    snakemake -s lf_Snakefile \
        --executor slurm \
        --default-resources slurm_partition=broadwl slurm_account=pi-yangili1 \
        --jobs 50 --use-conda --rerun-incomplete --retries 3
    cd ..
done
