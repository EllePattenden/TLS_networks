#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --job-name=clean_ref_min
#SBATCH --output=/fred/oz139/elle/log/floyds_mix_clean/%j
#SBATCH --ntasks=1 
#SBATCH --array=1-2000          
#SBATCH --mem-per-cpu=5000

# run from /fred/oz139/elle/    with 
# chunk=X  (used here to overcome limit in # of array jobs)
## then sbatch --job-name=clean_ref_min$chunk --export=chunk=$chunk process_floyd.sh

SCRATCH_DIRECTORY=/fred/oz139/elle/working/${SLURM_JOBID}
CODE_DIRECTORY=/fred/oz139/elle/code
OUTPUT_DIRECTORY=/fred/oz139/elle/output/floyds_mix/cleaned_sim_results
mkdir -p ${SCRATCH_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}
mkdir -p ${OUTPUT_DIRECTORY}/cleaned
mkdir -p ${OUTPUT_DIRECTORY}/cleaned/networks
mkdir -p ${OUTPUT_DIRECTORY}/to_smoosh                # smooshing = step 2, then can delete this folder! 
cp ${CODE_DIRECTORY}/min_candp_inp.R ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.1.1
echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID}"
Rscript min_candp_inp.R
cp ${SCRATCH_DIRECTORY}/summary*.RData ${OUTPUT_DIRECTORY}/cleaned
cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}/to_smoosh/
cp ${SCRATCH_DIRECTORY}/*.png  ${OUTPUT_DIRECTORY}/cleaned/networks/
rm -rf ${SCRATCH_DIRECTORY}