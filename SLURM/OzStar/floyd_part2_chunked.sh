#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --job-name=refute_min_actuallly
#SBATCH --output=/fred/oz139/elle/log/floyds_mix/%j
#SBATCH --ntasks=1 
#SBATCH --array=1-2000          
#SBATCH --mem-per-cpu=8000

# chunked version of code to get around max_array 
# run from /fred/oz139/elle/    with 
# chunk=4     # a number 1-33 ** number 33 only has 260 rows, so ammend when get around to that! 
## then sbatch --job-name=refute_min_$chunk --export=chunk=$chunk floyd_part2_chunked.sh

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}
CODE_DIRECTORY=/fred/oz139/elle/code/
OUTPUT_DIRECTORY=/fred/oz139/elle/output/floyds_mix/sim_results
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R ${SCRATCH_DIRECTORY}
cp /fred/oz139/elle/output/floyds_mix/test_combo_chunked.Rda ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
T="$(date +%s)"
echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID}"
Rscript part2_chunked.R
ls ${SCRATCH_DIRECTORY}/replicating_min*.Rdata
cp ${SCRATCH_DIRECTORY}/replicating_min*.Rdata  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "run # ${SLURM_ARRAY_TASK_ID} took ${T} seconds!"
rm -rf ${SCRATCH_DIRECTORY}

