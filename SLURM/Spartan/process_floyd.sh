#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --job-name=clean_ref_min
#SBATCH --output=/data/gpfs/projects/punim1783/log/floyds_mix_clean/%A_%a.out
#SBATCH --array=1-2000         
#SBATCH --mem-per-cpu=1000

# 25/01/23 code has been changed for actually replicating Min, need to change back if re-running with rho sep! 
# see comments below! 

# cd /data/gpfs/projects/punim1783/jobs
# version="actually_rep"        | actually_rep OR rho_sep 
# sbatch --export=version=$version process_floyd.sh 

for chunk in {1..22};   # (shouldn't need to do last 840 files!, but should get errors in array 2000)
do 
    SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${chunk}_${SLURM_ARRAY_TASK_ID}
    CODE_DIRECTORY=/data/gpfs/projects/punim1783/code
    if [ "$version" == "actually_rep" ]; then  
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/actually_replicating_min
    else 
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results
    fi
    mkdir -p ${SCRATCH_DIRECTORY}
    mkdir -p ${OUTPUT_DIRECTORY}
    mkdir -p ${OUTPUT_DIRECTORY}/cleaned
    mkdir -p ${OUTPUT_DIRECTORY}/to_smoosh                            # smooshing = step 2, then can delete this folder! 
    cp ${CODE_DIRECTORY}/process_min_inP.R ${SCRATCH_DIRECTORY}       
    cd ${SCRATCH_DIRECTORY}
    module load r/4.0.0 
    T="$(date +%s)"
    echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID} "
    Rscript process_min_inP.R ${chunk} ${version}                           
    cp ${SCRATCH_DIRECTORY}/summary*.RData ${OUTPUT_DIRECTORY}/cleaned/
    cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}/to_smoosh/
    T="$(($(date +%s)-T))"
    echo "run # ${SLURM_ARRAY_TASK_ID} in chunk ${chunk} took ${T} seconds"
    rm -rf ${SCRATCH_DIRECTORY}
done
