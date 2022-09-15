##--------------------------------------------------------------------------
##                      REPLICATE/REFUTE MIN ET AL.                       --
##                      Data Cleaning and Analysis                        --
##                               for OZSTAR                               --
##--------------------------------------------------------------------------

# GOALS: (1) recreate figure 3 from Min et al., while (2) getting rid of 
# the big data problem I accidentally created in the process. 
# ( at present, network objects are saved whenever a tie or strategy is 
#  updated; that's a lot of networks...)

# APPROACH: clean and process in parallel, then smoosh 'out' results; 
# delete network objects, unless they contain something *interesting* 

library(data.table)
library(tidyverse)
library(igraph)
library(Matrix)

rm(list = ls())           # clear environment 
options(scipen=999)       # stop printing in scientific notation, please !
igraph_options(sparsematrices = TRUE)  # work with sparse matrices, please

chunk <- Sys.getenv("chunk")              # processes batches of 2000 at a time
chunk <- chunk - 1                        # using 0 index 
k <- Sys.getenv("SLURM_ARRAY_TASK_ID")    # identifies file within ^ 

out <- data.table(rep = 0,             # output for identification and plotting 
                  ri = 0, 
                  rho = 0,
                  degree = 0,
                  f_prop_coop = 0,       # final proportion of coop.
                  ticked = 0,            # when the run ended
                  av_prop_coop = 0,      # average prop coop. over the run
                  sd_prop_coop = 0,      # sd of ^ 
                  num_tie_swaps = 0,
                  num_strat_swaps = 0
                  )
prop_ties <- data.table(cc = rep(NaN, 100001),    # "take-homes" from (boring) networks 
                        cd = rep(NaN, 100001),    # so that we can cut down
                        dd = rep(NaN, 100001))    # output size 

# ------------------------------ Get Data ------------------------------------#

files <- list.files("/fred/oz139/elle/output/floyds_mix/sim_results", 
                    recursive = FALSE,      # no sub folders 
                    pattern = ".Rdata",   
                    full.names = FALSE      # just store file names, not path
                    )

doing <- files[chunk + k]   # identify file to process 
data <- get(load(doing))    # load it into working environment 
                        # three components to 'data' list:  
output <- data[[1]]         # 1. summary stats for each tick
network_list <- data[[2]]   # 2. networks at t0, tf, and every tick a tie changed
debug <- network_list
node_list <- data[[3]]      # 3. attribute data 
rm(data, files)                    

# ------------------------------ Process -------------------------------------#

# Identify the run! 
# file names: "replicating_min_rep_XX_degree_XX_RI_XX_rho_XX_.Rdata"
# get those XXs and store in out 
ri <- str_match(doing, "RI_\\s*(.*?)\\s*_rho")[,2]
out[, 'ri'] <- ri
rho <- str_match(doing, "_rho_\\s*(.*?)\\s*_.R")[,2]
out[, 'rho'] <- rho
deg <- str_match(doing, "degree_\\s*(.*?)\\s*_RI")[,2]
out[, 'degree'] <- deg
rep <- str_match(doing, "rep_\\s*(.*?)\\s*_degree")[,2]
out[, 'rep'] <- rep

# fill in more of out 
output[, tick := .I - 1]                            # label time steps
output <- na.omit(output, "prop_coop")              # drop empties
out[, 'ticked'] <- output[which.max(tick), tick]    # store final tick # 
out[, 'f_prop_coop'] <- output[which.max(tick), prop_coop]  # store final p_coop
out[, 'av_prop_coop'] <- sum(output[-1,"prop_coop"]) / output[which.max(tick), tick]
out[, 'sd_prop_coop'] <- output[-1, sd(prop_coop)]
out[, 'num_tie_swaps'] <- sum(output[tie_updated == 1, tie_updated])
out[, 'num_strat_swaps'] <- sum(output[strat_updated == 1, strat_updated])

# process networks - 1. get breakdown for prop_ties 
network_list[sapply(network_list, is.null)] <- NULL         # remove empty
names(network_list) <- c(1, output[tie_updated == 1 | strat_updated == 1, tick]) # label
tally <- 1
for (i in network_list) {
  copy <- as.matrix(i)
  which <- names(network_list)[tally]
  strat <- node_list[tick == as.numeric(which), strategy]
  prop_ties[as.numeric(which), 'cc'] <- sum(copy[strat == 1, strat == 1])/2 
  prop_ties[as.numeric(which), 'cd'] <- sum(copy[strat == 1, strat == 0]) 
  prop_ties[as.numeric(which), 'dd'] <- sum(copy[strat == 0, strat == 0])/2
  tally <- tally + 1
  rm(copy, which, strat)
}
rm(tally)     
# fill in gaps 
copy<- as.matrix(network_list[[1]])  # starting with t0
strat<- node_list[tick == 0, strategy]
prop_ties[1, 'cc'] <- sum(copy[strat == 1, strat == 1])/2 
prop_ties[1, 'cd'] <- sum(copy[strat == 1, strat == 0]) 
prop_ties[1, 'dd'] <- sum(copy[strat == 0, strat == 0])/2
prop_ties[prop_ties == 'NaN' ] <- NA
total <- sum(prop_ties[1,])
prop_ties[, cc := cc[1], by = cumsum(!is.na(cc))][, cc:= cc / total]
prop_ties[, cd := cd[1], by = cumsum(!is.na(cd))][, cd := cd / total]
prop_ties[, dd := dd[1], by = cumsum(!is.na(dd))][, dd := dd / total]
prop_ties[, time := .I - 1]
# store summary plot of proportion over time 
prop_ties_time <- ggplot(prop_ties, aes(x = time)) + 
  geom_line(aes(y = cc), color = "green") + 
  geom_line(aes(y = cd), color = "orange") + 
  geom_line(aes(y = dd), color = "red") + 
  theme_classic() + ylab("proportion of ties") + 
  ggtitle(paste0("rep ", rep, " degree ", deg, " rho ", rho, " resource inflow ", ri))
prop_ties_time
prop_c_time <- ggplot(output, aes(x = tick)) +        # and for prop cooperating 
  geom_line(aes(y = prop_coop), color = "green") +     
  ylab("proportion cooperating") + theme_classic() + 
  ggtitle(paste0("rep ", rep, " degree ", deg, " rho ", rho, " resource inflow ", ri))
prop_c_time

# process networks - 2. create initial and final networks 
networks2plot <- c(head(names(network_list), 1), tail(names(network_list), 1))
for (net in networks2plot) {
  if (net == "1") {
    edge_list <- node_list[tick == (as.numeric(net) - 1), c("id", "strategy")]
  } else { 
    edge_list <- node_list[tick == (as.numeric(net)), c("id", "strategy")] 
  }
  network <- graph.adjacency(network_list[as.character(net)][[1]], mode = "undirected")
  V(network)$strategy <-edge_list$strategy
  V(network)$color <- ifelse(V(network)$strategy == 1, "green", "red")
  set.seed(3113)
  # plot and save
  png <- file.path(getwd(), 
                   paste0("rep_", rep, "_net_", net, ".png"))
  png(png)
  plot(network,
       layout = layout_with_fr,
       vertex.size = 10,
       vertex.label.cex = 0.5)
  dev.off()
  rm(png, network, edge_list)
}
# export and clean up 
name <- paste0(paste("cleanedNW","refMin", "degree", deg, "ri", ri, "rho", rho,"rep", rep, sep = "_"),".RData")
save(prop_ties, prop_ties_time, prop_c_time, file = name)
rm(name)
name <- paste0(paste("smoosh", "degree", deg, "ri", ri, "rho", rho,"rep", rep, sep = "_"), ".RData")
save(out, file = name)
rm(name)


