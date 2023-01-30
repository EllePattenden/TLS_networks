#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --job-name=clean_ref_min
#SBATCH --output=/data/gpfs/projects/punim1783/log/floyds_mix_clean/PROCESSING%A_%a.out
#SBATCH --array=1-2000      
#SBATCH --mem-per-cpu=1000

# cd /data/gpfs/projects/punim1783/jobs
# version="actually_rep"       | version="rho_sep"
# go="original"                | go="r50"
# sbatch --export=version=$version,go=$go process_floyd.sh 

CODE_DIRECTORY=/data/gpfs/projects/punim1783/code
if [[ ${version} == "actually_rep" ]]; then 
    OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/actually_replicating_min/
    mkdir -p ${OUTPUT_DIRECTORY}
    if [[ ${go} == "original" ]]; then 
        mkdir -p ${OUTPUT_DIRECTORY}/to_smoosh/ 
    else 
        # {go} == "r50"
        mkdir -p ${OUTPUT_DIRECTORY}/to_smoosh_R50/
    fi 
else 
    # ${version} == "rho_sep"
    OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results
    mkdir -p ${OUTPUT_DIRECTORY}
    if [[ ${go} == "original" ]]; then
        mkdir -p -v ${OUTPUT_DIRECTORY}/to_smoosh2/      
    else 
        # {go} == "r50"
        mkdir -p ${OUTPUT_DIRECTORY}/to_smoosh_R50/
    fi
fi
mkdir -p -v ${OUTPUT_DIRECTORY}/cleaned/

for chunk in {1..22};   # don't need to do last 840 files!; should get errors in array 2000
do 
    SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${chunk}_${SLURM_ARRAY_TASK_ID}
    mkdir -p -v ${SCRATCH_DIRECTORY}
    cp ${CODE_DIRECTORY}/process_min_inP.R ${SCRATCH_DIRECTORY}  
    cd ${SCRATCH_DIRECTORY}
    module load r/4.0.0 
    T="$(date +%s)"
    echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID} "
    Rscript process_min_inP.R ${chunk} ${version}   
    cp ${SCRATCH_DIRECTORY}/summary*.RData ${OUTPUT_DIRECTORY}/cleaned/                        
    if [[ ${go} == "r50" ]]; then 
        cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}/to_smoosh_R50/
    else 
        if [[ "${version}" == "actually_rep" ]]; then
            cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}/to_smoosh/
        else
            cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}/to_smoosh2/  
        fi
    fi
    T="$(($(date +%s)-T))"
    echo "run # ${SLURM_ARRAY_TASK_ID} in chunk ${chunk} took ${T} seconds"
    rm -rf ${SCRATCH_DIRECTORY}
done
