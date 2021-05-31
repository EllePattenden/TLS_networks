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
##  except when close to 1. However, Min et al. created a dependancy (i.e.,   
##  negative relationship) between the probability of strategy updating and   
##  probability of network rewiring, as set by the parameter rho. Here, we    
##  examine what happens when the probability of network rewiring is           
##  manipulated independently.   

# CONTEXT
##  This is part of a body of work that argues that, in conditions where   
##  collective action problems are truly large scale and/or involve the    
##  management of open access resources, a propensity to sever ties with   
##  defectors (which is advantageous for dyadic cooperation as well as     
##  scenarios where group-level selection dominates) is not be doing us   
##  any favours.     

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
options(scipen=999)  # stop printing in scientific notation !
set.seed(54321)      # need to swap if running in parallel ! 
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
                   # between cooperatorsand defectors with small resource 
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

# Systematically Varied 
# replicate <- seq(1, 30, 1)         # 30 reps of each
# degree <- c(15, 90)                # comparing average degree of 15 and 19
# resource_inflow <- seq(10, 60, 1)  # across c [10:60]
# rho <- sort(seq(0, 1, 0.01), decreasing = TRUE) # with rho [0:1]

#debug
replicate <- 1         
degree <- 15                
resource_inflow <- seq(10, 60, 10)  
rho <- sort(seq(0, 1, 0.1), decreasing = TRUE) # with rho [0:1]


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

modelRun(combo)

# ------------------------------ Get Data ------------------------------------#

files <- list.files(here("replicating_min"), 
                    recursive = TRUE,
                    pattern = ".Rdata",
                    full.names = TRUE)

data <- sapply(files, function(x) get(load(x)), simplify = FALSE)  
# now have a fat list, where each item contains the output (itself as a list)
list2env(data, envir = globalenv()) # stored seperately, but names shitttt 
rm(files, data)
for (i in ls()) {
  old <- as.list(i)
  new_name <- stringr::str_sub(old[[1]], start = 117)
  assign(sub(old[[1]], new_name, i), get(i))
  rm(i, old, new_name)
}
# this is not reproducible, but will work for now
rm(list = ls(pattern = '/Users/ellepattenden/Dropbox/PhD/2021/simulations/'))

# ------------------------------ Prelim ------------------------------------#

prop_coop_at_fin <- data.table(ri = rep(seq(10, 60, 10), times = 11), 
                               rho = rep(seq(0, 1, .1), each = 6),
                               prop_coop = NA,
                               ticked = NA)
prop_coop_at_fin[, names(prop_coop_at_fin) := lapply(.SD, as.numeric)]

stuff <- ls()[! ls() %in% "prop_coop_at_fin"]
for (i in stuff) { 
  print(i)
  # identify condition
  ri1 <- str_match(i, "RI_\\s*(.*?)\\s*_rho")
  ri1 <- ri1[, 2]
  rho1 <- str_match(i, "_rho_\\s*(.*?)\\s*_.Rdata")
  rho1 <- rho1[, 2]
  # get main output
  dat <- get(i)[[1]]
  dat[, tick := .I - 1]
  dat <- na.omit(dat, "prop_coop")
  # get ties
  ties <- get(i)[[2]]
  names(ties) <- seq_along(ties) - 1 # label with tick
  ties[sapply(ties, is.null)] <- NULL   # remove empty
  networks_to_plot <- names(ties)  # for plotting
  networks_to_plot <- c(head(networks_to_plot, 1),
                        tail(networks_to_plot, 1))
  # get attributes
  att <- get(i)[[3]]

  if (dat[which.max(prop_coop), prop_coop] != 1) {
    min_coop <- dat[which.min(prop_coop), prop_coop]
    ticked1 <- dat[which.min(prop_coop), tick]
    prop_coop_at_fin [ri == ri1 & rho == rho1,
                      prop_coop := min_coop]
    prop_coop_at_fin [ri == ri1 & rho == rho1, 
                      ticked := ticked1]
    rm(min_coop, ticked1)
  } else {
    max_coop <- dat[which.max(prop_coop), prop_coop]
    prop_coop_at_fin[ri == ri1 & rho == rho1,
                      prop_coop := max_coop]
    ticked1 <- dat[which.max(prop_coop), tick]
    prop_coop_at_fin [ri == ri1 & rho == rho1, 
                      ticked := ticked1]
    rm(max_coop, ticked1)
  }
  
  ggplot(dat, aes(x = tick)) +
    geom_line(aes(y = prop_coop)) +
    labs(title = as.character(i)) + ylab("proportion cooperating")

  for (net in networks_to_plot) {
    edge_list <- att %>% filter(tick == net) %>% select(id, strategy)
    network <- graph.adjacency(ties[as.character(net)][[1]],
                               mode = "undirected")
    V(network)$strategy <-  edge_list$strategy
    V(network)$color <- ifelse(V(network)$strategy == 1, "green", "red")
    set.seed(1234)
    plot.igraph(network, layout = layout_with_lgl,
                vertex.size = 10,vertex.label.cex = 0.5)
  }
  rm(dat, att, ties, ri1, rho1)
}
rm(stuff, i, edge_list, network, networks_to_plot)

# only working with limited number of runs at the moment
prop_coop_at_fin <- na.omit(prop_coop_at_fin, "ticked")

ggplot(prop_coop_at_fin, aes(x = ri, y = rho, fill = prop_coop)) + 
  geom_tile() + xlab("resource inflow, c") + ylab("rewiring probability, rho")