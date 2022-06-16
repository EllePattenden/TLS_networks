##--------------------------------------------------------------------------
##                      REPLICATE/REFUTE MIN ET AL.                       --
##--------------------------------------------------------------------------

# OVERVIEW
##  This set of simulations explore the purported relationship between tie       
##  rewiring probability and resoure inflow on a population of cooperating     
##  agents' resistance to invasion by a lone defector in an agent-based        
##  version of the TLS model (2012; 2016) with network dynamics described by   
##  Min et al. (2018). They found that increasing the probability of network   
##  rewiring had a positive effect in networks with low average dgree (k =     
##  15), but had minimal impact in relatively high degree networks (k = 90)       
##  except when close to 1. However, Min et al. created a dependency (i.e.,   
##  negative relationship) between the probability of strategy updating and   
##  probability of network rewiring, as set by the parameter rho. Here, we    
##  examine what happens when the probability of network rewiring is           
##  manipulated independently.

## A year later - Elle realises that we will have to do the follow ups, (1) hold  
# rho constant and vary the probability of strategy updating and (2) vary both as
# the the prob of strategy updating is tied to the payoff difference, which is 
# tied to the resource inflow... 

# CONTEXT
##  This is part of a body of work that argues that, in conditions where   
##  collective action problems are truly large scale and/or involve the    
##  management of open access resources, a propensity to sever ties with   
##  defectors (which is advantageous for dyadic cooperation as well as     
##  scenarios where group-level selection dominates) may not be doing us   
##  any favours.     

# https://dashboard.hpc.unimelb.edu.au/forms/script_generator/

##-------------------------------------------------------------------------#

# install.packages(c("data.table", "tidyverse", "igraph", "rlist",
#                    "Matrix", "here"))
library(data.table)
library(tidyverse)
library(igraph)
library(rlist)
library(Matrix)
library(here)

here::here()         # set directory
rm(list = ls())      # clear environment 
options(scipen=999)  # stop printing in scientific notation, please !
set.seed(31133113)      # swap so each make brain do sense thing and batch (elle think)
igraph_options(sparsematrices = TRUE)  # work with sparse matrices, please

# -------------------------- Parameters ----------------------------------#

# Fixed 
params <- setDT(data.frame(
  N = 100,   # number of agents 
  # Resource related
  initi_R = 50,   # initial size of resource
  max_R = 200,    # resource carrying capacity
  d = 50,         # natural discharge rate
  # Resource 
  E_opt = 0.483,  # optimal effort (water)
  q = 1,          # technology quotient 
  w = 15,         # opportunity cost per unit effort 
  # Production (Cobb-Douglas)
  alpha = 0.6,    
  beta = 0.2,     # with a + b < 1 
  gamma = 10,     # for diminishing returns to scale 
  # Ostracism (Gompertz function)
  h = .17,        # maximum sanctioning (asymptote)      
                  # Noting that all papers bar Min et al. use ~.333 !!! 
  t = -150,       # sanctioning effectiveness threshold
  g = -10,        # sanctioning growth rate
  # Other 
  max_t = 10^5,    # max time steps in each  simulation i.e., time for
                   # the system to equilibrate; 50,000 in TLS et al. (2016)
                   # 100,000 in Min et al. (2018) because ties continue to be
                   # updated until homogeneous communities form 
                   # will break this is subsequent sims by using utility comps.
  p_coop = 0.99,   # initial proportion of cooperators
  mu = 3.9,        # "cheating multiplier"; 0 < mu < 3.9, where 3.9 = Nash
                   # (overwritten for a more precise value below) 
  negative_utilities = TRUE,   # can agent utilities be negative? logical;
                   # false = they can't. TLS et al. (2016) said no 
                   # but you * might * need  them here to get differences  
                   # between cooperator sand defectors with small resource 
                   # inflows here... (Min et al. allowed, I think)
  update_prob = 0.01,  # proportion of agents comparing utilities in each tick 
                   # referred to as "strategy updating probability" in some
  network_type = "ER_random",   # the initial network structure
  scheduling_dynamics = "seperate",  # how tie formation plays out. 
                  # For 'seperate' (here), the probability of rewiring 
                  # probability is determined by rho while the probability 
                  # for strategy updating is proportional to the payoff 
                  # difference. 
  tie_strat = "min_method"  # how agent's rewire their networks
                  # as this simulation seperates the prob. of strategy 
                  # updating from the prob of network rewiring. Two agents 
                  # are randomly selected each tick; one for each social 
                  # process.
  # noting that the last 5 parameters are set here for reference only. 
))
params[, "e_coop" := E_opt / N ]    # effort for cooperators 
params[, "e_defect" := e_coop *  mu]    # defector's effort
# params[, "e_defect" := 1.826 / N]   # ' ' if you want the rounding to be exact 

# Systematically Varied = 309060 runs 
# replicate <- seq(1, 30, 1)         # 30 reps of each
# degree <- c(15, 90)                # comparing average degree of 15 and 19
# resource_inflow <- seq(10, 60, 1)  # across c [10:60]
# rho <- sort(seq(0, 1, 0.01), decreasing = TRUE) # with rho [0:1]

#debug / pilot 
# replicate <- 1
# degree <- 15
# resource_inflow <- seq(10, 60, 10)
# rho <- sort(seq(0, 1, 0.1), decreasing = TRUE) # with rho [0:1]

# compromise for Floyd = 64260 runs 
replicate <- seq(1, 30, 1)         # 25 reps of each
degree <- c(15, 90)                # comparing average degree of 15 and 90
resource_inflow <- seq(10, 60, 1)  # across c [10:60]
rho <- sort(seq(0, 1, 0.05), decreasing = TRUE) # with rho [0:1]

# stick all parameters in one big data.table 
combos <- setDT(expand.grid(replicate, degree, resource_inflow, rho))
setnames(combos, new = c("replicate", "degree", "resource_inflow", "rho"))
combo <- cbind(combos, params)   # each row has parameter values for a run 

# keep it clean(-ish) 
rm(list = 
     c("replicate", "degree", "resource_inflow", "rho", "combos", "params"))

# ----------------------- Set-Up Procedures --------------------------------#

# get functions to generate a population of agents:
source("agent_generator.R")   # agent_generator() and network_generator()

# get functions that run the simulations:            
source("one_run_min.R")   # one_run() - implements the go procedure
                          # modelRun() - runs one_run for all rows of combo

# set up folder for storing output 
if (!dir.exists("replicating_min")) {     
  dir.create("replicating_min")         # create main folder
} 

# ---------------------------- Run and Save ----------------------------------#

modelRun(combo)  # run for full parameter space

purrr::pmap(combo[1,], one_run)   # run for one set of parameter combinations



