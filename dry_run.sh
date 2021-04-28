#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH -J dryRun
#SBATCH -o dryRun.%N.%j.out
#SBATCH -e dryRun.%N.%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
####SBATCH -p unlimitq

# Environment
module purge
module load bioinfo/snakemake-5.8.1
module load system/singularity-3.6.4

# Variables
CONFIG=config/genologin.yaml
COMMAND="sbatch --nodelist={cluster.nodelist} --mem={cluster.mem} --ntasks-per-node {cluster.npernode} -t {cluster.time} -n {cluster.ntasks} -c {cluster.c} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err"
CORES=32
mkdir -p snake_subjob_log

# Workflow
snakemake -s Snakefile -n -r --use-singularity -j $CORES --cluster-config $CONFIG --cluster "$COMMAND"
