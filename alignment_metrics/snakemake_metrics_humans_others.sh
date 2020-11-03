#!/bin/bash
# Job name:
#SBATCH --job-name=human_others
#
# Account:
#SBATCH --account=fc_genomicdata
#
# Partition:
#SBATCH --partition=savio,savio2,savio3,savio_bigmem,savio2_bigmem
#
# Quality of Service:
#SBATCH --qos=savio_normal
#
# Requeue if cancelled
#SBATCH --requeue
#
# Wall clock limit:
#SBATCH --time=3-00:00:00
#
# Command(s) to run:
snakemake -s Snakefile_vst_variance_pbs_dcgh --configfile human_metrics_config.json --unlock
snakemake -s Snakefile_vst_variance_pbs_dcgh --configfile human_metrics_config.json --cluster "sbatch {params.slurm_opts}"  -j 6 --ri -k -w 120

















