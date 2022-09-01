#!/bin/bash
#SBATCH --time=00:05:00
#SBATCH --job-name=refute_min_test
#SBATCH --output=/fred/oz139/elle/log/floyds_mix/%j
#SBATCH --mem-per-cpu=5000
#SBATCH --cpus-per-task=2

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}/
CODE_DIRECTORY=/fred/oz139/elle/code
OUTPUT_DIRECTORY=/fred/oz139/elle/output/floyds_mix
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
cp ${CODE_DIRECTORY}/floyds_mix.R ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
Rscript floyds_mix.R
ls ${SCRATCH_DIRECTORY}
cp ${SCRATCH_DIRECTORY}/test_combo.Rda  ${OUTPUT_DIRECTORY}
rm -rf ${SCRATCH_DIRECTORY}