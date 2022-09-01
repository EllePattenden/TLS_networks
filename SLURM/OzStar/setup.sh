#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --job-name=refute_min_test
#SBATCH --output=/fred/oz139/elle/log/%j
#SBATCH --mem-per-cpu=5000
#SBATCH --cpus-per-task=2

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}/
CODE_DIRECTORY=/fred/oz139/elle/code/
OUTPUT_DIRECTORY=/fred/oz139/elle/output/
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
T="$(date +%s)"
Rscript testing.R
ls ${SCRATCH_DIRECTORY}
cp ${SCRATCH_DIRECTORY}/*.Rd*  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "It took ${T} seconds!"
rm -rf ${SCRATCH_DIRECTORY}