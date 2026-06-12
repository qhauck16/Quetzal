#!/usr/bin/env bash
# EXAMPLE per-chromosome dispatcher for the leafcutter -> fastTopics step.
#
# v0.1 deliberately ships NO required dispatch wrapper -- parallelise
# across chromosomes however your cluster prefers. This file is one
# concrete pattern that worked for the manuscript: for every chr<N>/
# subdir present under data/all_genes/, materialise a sibling chr<N>/
# directory next to lf_Snakefile, copy the Snakefile in, and submit
# one Snakemake job (SLURM, --use-conda) restricted to that chromosome.
#
# Fill in the two SLURM placeholders below (or rewrite the snakemake
# invocation entirely for your scheduler). v1.0 will provide a built-in
# dispatcher; for now this is just a starting point you can copy.
#
# Run from <repo>/scripts/genome_wide/:
#     bash setting_up_snakemake.sh

set -euo pipefail

ALL_GENES_DIR="../../data/all_genes"

# >>>>> EDIT FOR YOUR CLUSTER <<<<<
SLURM_PARTITION="<your-slurm-partition>"
SLURM_ACCOUNT="<your-slurm-account>"
MAX_JOBS=50

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
        --default-resources slurm_partition="$SLURM_PARTITION" slurm_account="$SLURM_ACCOUNT" \
        --jobs "$MAX_JOBS" --use-conda --rerun-incomplete --retries 3
    cd ..
done
