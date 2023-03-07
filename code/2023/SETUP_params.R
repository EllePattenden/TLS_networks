##--------------------------------------------------------------------------
##                           SETUP parameters                             --
##                      TLS with DYNAMIC NETWORKS                         --
##--------------------------------------------------------------------------

# It's 2023, time to get our shit together ...

# -------------------------- Parameters -----------------------------------#

# install.packages(c("data.table","tidyverse","igraph","rlist","Matrix",))

library(data.table)
library(tidyverse)
library(igraph)
library(Matrix)

rm(list = ls())       # clear environment 
options(scipen=999)   # stop printing in scientific notation, please !
igraph_options(sparsematrices = TRUE)  # work with sparse matrices
set.seed(31133113)   

# Fixed 
params <- setDT(data.frame(
  N = 100,   # number of agents 
  # Resource related
  initi_R = 50,   # initial size of resource
  max_R = 200,    # resource carrying capacity
  d = 50,         # natural discharge rate
  E_opt = 0.483,  # optimal effort (water use)
  q = 1,          # technology quotient 
  w = 15,         # opportunity cost per unit effort 
  # Production (Cobb-Douglas)
  alpha = 0.6,    
  beta = 0.2,     # with a + b < 1 
  gamma = 10,     # for diminishing returns to scale 
  # Ostracism (Gompertz function)
  h = .17,        # maximum sanctioning (asymptote)      
                  # Noting that all papers bar Min et al. 
                  # use ~.333 !!! 
  t = -150,       # sanctioning effectiveness threshold
  g = -10,        # sanctioning growth rate
  # Other 
  max_t = 10^5,   # max time steps in each  simulation 
                  # i.e., time for the system to equilibrate;
                  # 50,000 in TLS et al. (2016)
                  # 100,000 in Min et al. (2018) 
  p_coop = 0.99,  # initial proportion of cooperators
  mu = 3.9,       # "cheating multiplier"; 0 < mu < 3.9, 
                  # where 3.9 = Nash
                  # (can be overwritten for a more precise value below) 
  negative_utilities = TRUE,   # Can agent utilities be negative? Logical:
                               # F = they can't. TLS et al. (2016) said no 
                               # (Min et al. allowed, I think)
  update_prob = 0.01,  # proportion of agents comparing utilities in each tick 
                       # NOTE, this isn't the "strategy updating probability" 
  network_type = "ER_random"   # the initial network structure
))

params[, "e_coop" := E_opt / N ]       # effort for cooperators 
params[, "e_defect" := e_coop *  mu]   # defector's effort
# params[, "e_defect" := 1.826 / N]    # ' ' if you want the rounding to be exact 

# parameters updated later on 
params[, "version" := "hold"] 
params[, "network_dynamics" := "hold"]
params[, "ostracism_type" := "hold"]

# Parameter space, compromising for Floyd 
# = 42840 runs 
replicate <- seq(1, 20, 1)         # 20 reps of each
degree <- c(15, 90)                # comparing average degree of 15 and 90
resource_inflow <- seq(10, 60, 1)  # across c [10:60]
rho <- sort(seq(0, 1, 0.05), decreasing = TRUE) # with rho [0:1]
  
# stick all parameters in one big data.table and save as parameters.Rda
combos <- setDT(expand.grid(replicate, degree, resource_inflow, rho))
setnames(combos, 
         new = c("replicate", "degree", "resource_inflow", "rho"))
combo <- cbind(combos, params)     
save(combo, file="parameters.Rda") # each row has parameter values for a run

# create chunked version and save as test_combo_chunked.Rda
length <- ceiling(nrow(combo)/2000)   # chunks required 
out <- vector("list", length = length)   
for (i in 1:length(out)) { 
  start <- (i*2000) - 1999
  row <- i*2000
  out[[i]] <- combo[start:row, ]
}
# drop NA rows from chunk 22 
out[[22]] <- na.omit(out[[22]])
combo <- out
save(combo, file = "parameters_chunked.Rda")

# keep it clean(-ish) 
rm(list = c("i", "out", "row", "start", "length", "replicate", 
            "degree", "resource_inflow", "rho", "combos", "params"))

