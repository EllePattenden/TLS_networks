#!/bin/bash
#SBATCH --time=7:30:00
#SBATCH --job-name=refute_min_chunk
#SBATCH --output=/data/gpfs/projects/punim1783/log/floyds_mix/chunk20/%j  
#SBATCH --ntasks=1 
#SBATCH --array=1-2000                          
#SBATCH --mem-per-cpu=3500        # 2500 more than enough for low rho      

# cd /data/gpfs/projects/punim1783/jobs    
# chunk=XX                                  # a number 1-22 (noting chunk 22 only has 840 rows)
# version="actually_rep"      | version="rho_sep"  
# go="original"               | go="r50"
# echo "doing chunk: " $chunk " version: " $version "go procedure: " $go
# sbatch --job-name=refute_min_$chunk --export=chunk=$chunk,version=$version,go=$go floyd_part2_chunked.sh

SCRATCH_DIRECTORY=/data/gpfs/projects/punim1783/working/${SLURM_JOBID}
CODE_DIRECTORY=/data/gpfs/projects/punim1783/code/
if [[ ${version} == "actually_rep" ]]; then 
    if [[ ${go} == "original" ]]; then 
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results/toreproducemin/
    else 
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results/toreproducemin/R50/
    fi 
else 
    if [[ ${go} == "original" ]]; then
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results/
    else 
        OUTPUT_DIRECTORY=/data/gpfs/projects/punim1783/output/floyds_mix/sim_results/R50/
    fi
fi
mkdir -p -v ${OUTPUT_DIRECTORY}
mkdir -p -v /data/gpfs/projects/punim1783/log/floyds_mix/chunk${chunk}/
mkdir -p -v ${SCRATCH_DIRECTORY}
cp ${CODE_DIRECTORY}/*.R ${SCRATCH_DIRECTORY}
cp /data/gpfs/projects/punim1783/output/floyds_mix/test_combo_chunked.Rda ${SCRATCH_DIRECTORY}
cd ${SCRATCH_DIRECTORY}
module load r/4.0.0 
T="$(date +%s)"
echo "beginning chunk ${chunk} run ${SLURM_ARRAY_TASK_ID}"
Rscript part2_chunked.R  ${chunk} ${version} ${go}          # this script will use version and go to call the desired GO script 
ls ${SCRATCH_DIRECTORY}/*replicating_min*.Rdata
cp ${SCRATCH_DIRECTORY}/*replicating_min*.Rdata  ${OUTPUT_DIRECTORY}
T="$(($(date +%s)-T))"
echo "run # ${SLURM_ARRAY_TASK_ID} took ${T} seconds"
rm -rf ${SCRATCH_DIRECTORY}
