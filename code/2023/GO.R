##--------------------------------------------------------------------------
##                                 "GO"                                   --
##                      TLS with DYNAMIC NETWORKS                         --
##                         (less garbo version)                           --
##--------------------------------------------------------------------------

# It's 2023, let's get our shit together! 

# Script replaces part2_chunked.R (and one_run_min*.R), which was called by 
# floyd_part2_chunked.sh 

# Called by sims_chunked.sh and sources agent_generator.R , 
# ostracism_functions.R , network_dynamics_functions.R , and SETUP_GO.R 

##-------------------------------------------------------------------------#
# install.packages(c("data.table","tidyverse","igraph","Matrix",))
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(Matrix))
rm(list = ls())        # clear environment 
options(scipen=999)    # stop printing in scientific notation, please !     
igraph_options(sparsematrices = TRUE) 

# get arguments from script call in slurm script 
args <- commandArgs(trailingOnly=TRUE)    
chunk <- as.numeric(args[1])    # ${chunk}= chunk %in% 1:22
ver <- args[2]                  # ${version}= "actually_rep" | "rho_sep"
network_d <- args[3]            # ${network_dynamics}= "min_method" |"NIM"
ost_type <- args[4]             # ${ostracism_type}= "original" | "FEJ" 
k <- as.numeric(                # k %in% 1:2000 (row in chunk)
  Sys.getenv("SLURM_ARRAY_TASK_ID")) 
set.seed(k)                     # make results reproducible

load("parameters_chunked.Rda")  # reads in a list of data.tables, 
                                # from which the right chunk ^ and  
                                # row (k) needs to be extracted
combo <- combo[[chunk]] 
combo[, "version" := ver]      
combo[, "network_dynamics" := network_d]
combo[, "ostracism_type" := ost_type]
run <- as.data.frame(combo)

# source functions used in GO! 
source("agent_generator.R")  
source("ostracism_functions.R") 
source("network_dynamics_functions.R")

if(ver == "rho_sep") {
  source("SETUP_GO.R")            
} else {   # version == "actually_rep
  # I haven't compartmentalized these scripts because I don't 
  # think I'll be using them again... 
  if (network_d == "original") {
    source("/data/gpfs/projects/punim1783/code/one_run_min_actual.R")
  }
  if (network_d == "NIM")  {
    print("Elle, dummy, you never did this...")
  }
}
rm(args, ver, network_d, ost_type)

# Actually GO! 
purrr::pmap(run[k,], one_run)
