#!/usr/bin/env sh

######################################################################
# @author      : qhauck (qhauck@midway2-login2.rcc.local)
# @file        : setting_up_snakemake
# @created     : Thursday Jan 23, 2025 14:47:54 CST
#
# @description : 
######################################################################

for file in ../../all_genes/chr*.tsv
do
    chr=$(basename $file .tsv )
    mkdir -p $chr
    
    cp lf_Snakefile $chr/

    cd $chr
    snakemake -s lf_Snakefile --executor slurm --default-resources slurm_partition=broadwl slurm_account=pi-yangili1 --jobs 50 --use-conda --rerun-incomplete --retries 3 
    cd ..
done

