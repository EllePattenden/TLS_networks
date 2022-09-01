#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --job-name=refute_min_chunk
#SBATCH --output=/data/gpfs/projects/punim1783/log/floyds_mix/%j
#SBATCH --ntasks=1 
#SBATCH --array=1-2000          
#SBATCH --mem-per-cpu=4000            # TEST MEMORY REQS NOW THAT WE AREN'T SAVING ALL THE BLOODY NETWORKS

# chunked version of code to get around max_array 
# cd /data/gpfs/projects/punim1783/jobs    
# chunk=1     # a number 1-22 ** number 22 only has 840 rows, so ammend when get around to that! 
## then sbatch --job-name=refute_min_$chunk --export=chunk=$chunk floyd_part2_chunked.sh


SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${SLURM_JOBID}
CODE_DIRECTORY=/data/gpfs/projects/punim1783/code
OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R ${SCRATCH_DIRECTORY}
cp /data/gpfs/projects/punim1783/output/floyds_mix/test_combo_chunked.Rda ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.0.0 
T="$(date +%s)"
echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID}"
Rscript part2_chunked.R
ls ${SCRATCH_DIRECTORY}/replicating_min*.Rdata
cp ${SCRATCH_DIRECTORY}/replicating_min*.Rdata  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "run # ${SLURM_ARRAY_TASK_ID} took ${T} seconds"
rm -rf ${SCRATCH_DIRECTORY}
