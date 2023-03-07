#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --job-name=clean_ref_min
#SBATCH --output=/data/gpfs/projects/punim1783/log/2023/processing/PROCESSING_%A_%a.out
#SBATCH --array=1-2000      
#SBATCH --mem-per-cpu=1000

# version="rho_sep"              "rho_sep" | "actually_rep"
# network_dynamics="NIM"         "min_method" | "NIM"
# ostracism_type="original"      "original" | "FEJ"
# sbatch --export=version=$version,network_dynamics=$network_dynamics,ostracism_type=$ostracism_type /data/gpfs/projects/punim1783/jobs/2023/process_floyd.sh 

# version="rho_sep" network_dynamics="NIM" ostracism_type="original"

CODE_DIRECTORY=/data/gpfs/projects/punim1783/code/2023
OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/2023/cleaned_output/${version}/${network_dynamics}/${ostracism_type}
RAW_DIRECTORY=/data/gpfs/projects/punim1783/output/2023/raw_output/${version}/${network_dynamics}/${ostracism_type}

# to. calc extra stuff from previous runs ! 
# version="rho_sep" network_dynamics="min_method" ostracism_type="original" 
# RAW_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results
# version="actually_rep" network_dynamics="min_method" ostracism_type="original" 
# RAW_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results/toreproducemin


mkdir -pv ${OUTPUT_DIRECTORY}

for chunk in {1..22};   # don't need to do last 840 files!; should get errors in array 2000
do 
    SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${chunk}_${SLURM_ARRAY_TASK_ID}
    mkdir -pv ${SCRATCH_DIRECTORY}
    cp ${CODE_DIRECTORY}/process_min_inP.R ${SCRATCH_DIRECTORY}  
    cd ${SCRATCH_DIRECTORY}
    module load r/4.0.0 
    T="$(date +%s)"
    echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID} "
    Rscript process_min_inP.R ${chunk} ${version} ${network_dynamics} ${ostracism_type} ${RAW_DIRECTORY}
    cp ${SCRATCH_DIRECTORY}/smoosh*.RData ${OUTPUT_DIRECTORY}
    ls ${SCRATCH_DIRECTORY}/smoosh*.RData
    T="$(($(date +%s)-T))"
    echo "run # ${SLURM_ARRAY_TASK_ID} in chunk ${chunk} took ${T} seconds"
    rm -rf ${SCRATCH_DIRECTORY}
done
