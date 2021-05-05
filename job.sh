#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH -J detMut
#SBATCH -o detMut.%N.%j.out
#SBATCH -e detMut.%N.%j.err
#SBATCH --mem=1G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
####SBATCH -p unlimitq

# Environment
module purge
module load bioinfo/snakemake-5.25.0
module load system/singularity-3.6.4

# Variables
CONFIG=config/genologin.yaml
COMMAND="sbatch --cpus-per-task={cluster.cpus} --time={cluster.time} --mem={cluster.mem} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err"
# COMMAND="sbatch -p unlimitq --cpus-per-task={cluster.cpus} --time={cluster.time} --mem={cluster.mem} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err"
CORES=32
mkdir -p snake_subjob_log

# Workflow
snakemake -s Snakefile --use-singularity -j $CORES --cluster-config $CONFIG --cluster "$COMMAND"

## Session informations
echo '########################################'
echo 'Date:' $(date --iso-8601=seconds)
echo 'User:' $USER
echo 'Host:' $HOSTNAME
echo 'Job Name:' $SLURM_JOB_NAME
echo 'Job ID:' $SLURM_JOB_ID
echo 'Number of nodes assigned to job:' $SLURM_JOB_NUM_NODES
echo 'Nodes assigned to job:' $SLURM_JOB_NODELIST
echo 'Directory:' $(pwd)
## Detail Information:
#echo 'scontrol show job:'
#scontrol show job $SLURM_JOB_ID
echo '########################################'
