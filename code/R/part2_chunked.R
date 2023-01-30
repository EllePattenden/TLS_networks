# chunked = procedure for doing the dumb # of runs required to rep Min 

##--------------------------------------------------------------------------------------#
library(data.table)
library(tidyverse)
library(igraph)
library(Matrix)

rm(list = ls())        # clear environment 
options(scipen=999)    # stop printing in scientific notation, please !     
igraph_options(sparsematrices = TRUE) 

args <- commandArgs(trailingOnly=TRUE)   # ${chunk} ${version} ${go} 
chunk <- as.numeric(args[1])             # identifying chunk (1:22)
#chunk <- as.numeric(Sys.getenv("chunk")) 
version <- args[2]                       # version="actually_rep" | version="rho_sep"
go <- args[3]                            # go="original" | go="r50"

load("test_combo_chunked.Rda")     # now reads in a list of data.tables, from which  
                                   # the right chunk and then right row (k) needs 
                                   # to be extracted
source("agent_generator.R") 

    if(version == "actually_rep") {
        # for the actually replicating min runs
        if (go == "original") {
            source("one_run_min_actual.R")       
        } else {
            source("one_run_min_actual_R50.R")    # with schedule for when resource is updated shifted 
        }
    } else {   
        # version =="rho_sep" (i.e., refuting Min)
        if (go == "original") {
            source("one_run_min.R")
        } else {
            source("one_run_min_R50.R")          # with schedule for when resource is updated shifted 
        }
    }

k <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))    
set.seed(k) 

combo <- combo[[chunk]] 
run <- as.data.frame(combo)
purrr::pmap(run[k,], one_run)
