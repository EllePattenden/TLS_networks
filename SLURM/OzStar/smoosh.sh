#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=smoosh
#SBATCH --output=/fred/oz139/elle/log/floyds_mix_clean/%j
#SBATCH --ntasks=1     
#SBATCH --mem-per-cpu=5000

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}
CODE_DIRECTORY=/fred/oz139/elle/code
OUTPUT_DIRECTORY=/fred/oz139/elle/output/floyds_mix/cleaned_sim_results
mkdir -p ${SCRATCH_DIRECTORY}
cp ${CODE_DIRECTORY}/smoosh.R ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
Rscript smoosh.R
cp ${SCRATCH_DIRECTORY}/summary_min.RData ${OUTPUT_DIRECTORY}
rm -rf ${SCRATCH_DIRECTORY}
