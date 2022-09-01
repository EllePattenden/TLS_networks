#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --job-name=refute_min_actuallly
#SBATCH --output=/fred/oz139/elle/log/floyds_mix/%j
#SBATCH --ntasks=1 
#SBATCH --array=2047-3000        # 64260     -- LOOKS LIKE WE NEED TO DIVIDE COMBO UP! 
#SBATCH --mem-per-cpu=8000

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}
CODE_DIRECTORY=/fred/oz139/elle/code/
OUTPUT_DIRECTORY=/fred/oz139/elle/output/floyds_mix/sim_results
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R ${SCRATCH_DIRECTORY}
cp /fred/oz139/elle/output/floyds_mix/test_combo.Rda ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
T="$(date +%s)"
Rscript part2.R
ls ${SCRATCH_DIRECTORY}/replicating_min*.Rdata
cp ${SCRATCH_DIRECTORY}/replicating_min*.Rdata  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "run # ${SLURM_ARRAY_TASK_ID} took ${T} seconds!"
rm -rf ${SCRATCH_DIRECTORY}

