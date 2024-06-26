##--------------------------------------------------------------------------
##                      REPLICATE/REFUTE MIN ET AL.                       --
##                      Data Cleaning and Analysis                        --
##                              for SPARTAN                               --
##--------------------------------------------------------------------------

# GOALS: (1) recreate figure 3 from Min et al.(2018) and (2) identify runs 
# with *interesting* dynamics - relatively stable cooperation, followed by a 
# tipping point - that can be probed further. 

# APPROACH: clean and process in parallel, then smoosh 'out' results. 

##------------------------------- Setup ------------------------------------

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(igraph))
suppressMessages(library(Matrix))

rm(list = ls())                        # clear environment 
options(scipen=999)                    # go away scientific notation 
igraph_options(sparsematrices = TRUE)  # work with sparse matrices


args <- commandArgs(trailingOnly=TRUE)  # 2000 jobs run in p, each processing
# print(args)                             # 22 runs of the sim sequentially 
chunk <- as.numeric(args[1])            # identifying chunk (1:22)
version2 <- args[2]
network_dynamics <- args[3]
ostracism_type <- args[4]
directory <- args[5]
k <- as.numeric(                        # identify k (1:2000) 
  Sys.getenv("SLURM_ARRAY_TASK_ID"))     

out <- data.table(rep = 0,               # output for ID + plotting 
                  ri = 0,                # resource inflow
                  rho = 0,               # prob. strategy updating
                  degree = 0,            # average k: 15 or 90 
                  f_prop_coop = 0,       # final proportion of cooperators
                  ticked = 0,            # when the run ended   
                                           # 26/01/23 - now accurate for runs 
                                           # where stoppping rule was applied!  
                  av_prop_coop = 0,      # mean prop coop. over the run
                  sd_prop_coop = 0,      # sd of ^ 
                  num_tie_swaps = 0,     # tally times ties swapped
                  num_strat_swaps = 0,   # tally times strategies swapped,
                  final_CC = 0,          # proportion C-C ties at t_fin
                  final_CD = 0,          # proportion C-D ties at t_fin
                  final_DD = 0,          # proportion D-D ties at t_fin
                  num_strat_path_entered = 69696969,    # 26/01/23 tally times strat swap path entered - only relevant for the actually_replicating min runs!
                  final_water = 0,       # level of resource in last tick 
                  av_water = 0,          # average size of resource across run 
                  final_harvest = 0,     # size of total harvest in last tick 
                  av_harvest = 0,        # average harvest across run 
                  final_pi = 0,          # payoff (average) in final tick 
                  av_pi = 0,             # average payoff over run
                  final_utility = 0,     # utility (average)in final tick
                  av_utility = 0         # average utility over run
                  )

prop_ties <- data.table(cc = rep(NaN, 100001),    # to store summary
                        cd = rep(NaN, 100001),    # of networks over 
                        dd = rep(NaN, 100001))    # time 


# ------------------------------ Get Data ------------------------------------#

# should really be doing this in the slurm script but whatever... 
files <- list.files(directory, 
    recursive = FALSE,      # no sub folders 
    pattern = ".Rdata",   
    full.names = TRUE       # need the file path to load
)

# identify file to process 
doing <- files[ ((k - 1) * 22 + chunk )]       # (using 0 index)

# load it into working environment 
data <- get(load(doing))    
# there are 3 components to 'data' list:  
output <- data[[1]]         # 1. summary stats for each tick
network_list <- data[[2]]   # 2. t0 network + breakdown of C-C, C-D, D-D ties 
                            #    thereafter (with gaps, filled in below)
node_list <- data[[3]]      # 3. attribute data  (only used briefly) 
rm(data, files)                    

# ------------------------------ Process -------------------------------------#

# Identify the run! 
# file names: "XX_replicating_min_rep_XX_degree_XX_RI_XX_rho_XX_.Rdata"
# get those XXs and store in out 
ri <- str_match(doing, "RI_\\s*(.*?)\\s*_rho")[,2]
out[, 'ri'] <- ri
rho <- str_match(doing, "_rho_\\s*(.*?)\\s*_.Rdat")[,2]
out[, 'rho'] <- rho
deg <- str_match(doing, "degree_\\s*(.*?)\\s*_RI")[,2]
out[, 'degree'] <- deg
rep <- str_match(doing, "rep_\\s*(.*?)\\s*_degree")[,2]
out[, 'rep'] <- rep
print(                                                # read out for log
  paste("degree", deg, "| ri", ri, "| rho", rho,"| rep", rep, sep = " "))

# fill in more of out 
output[, tick := .I - 1]                               # label time steps
output <- na.omit(output, "prop_coop")                 # drop empties
out[, 'ticked'] <- output[which.max(tick), tick]       # store final tick # 
out[, 'f_prop_coop'] <- output[which.max(tick), 
                               prop_coop]              # final p_coop
out[, 'av_prop_coop'] <- sum(output[-1,"prop_coop"]) / 
                        output[which.max(tick), tick]  # av_p_coop
out[, 'sd_prop_coop'] <- output[-1, sd(prop_coop)]     # sd_p_coop
out[, 'num_tie_swaps'] <- sum(output[tie_updated == 1, # count tie swaps
                                     tie_updated])
out[, 'num_strat_swaps'] <- sum(output[strat_updated == 1, 
                                       strat_updated]) # count stat swaps
out[,'final_water'] <- output[which.max(tick), R_now]  # level of resource in last tick 
out[,'av_water'] <- sum(output[-1,R_now]) / 
                    output[which.max(tick), tick]      # av size of resource 
out[,'final_harvest'] <- output[which.max(tick), total_harvest]
out[,'av_harvest']   <- sum(output[-1,total_harvest]) / 
                        output[which.max(tick), tick]

# now dealing with NAs ...
out[,'final_pi'] <- output[which.max(tick), 
                      fifelse(!is.na(coop_U) & !is.na(defector_pi), (coop_U + defector_pi) /2, 
                      fifelse(is.na(coop_U), defector_pi, coop_U)) ]
out[,'av_pi'] <- sum(output[-1, fifelse(!is.na(coop_U) & !is.na(defector_pi), (coop_U + defector_pi),
                                fifelse(is.na(coop_U), defector_pi, coop_U))], 
                      na.rm=TRUE) / output[which.max(tick), tick]
out[,'final_utility'] <- output[which.max(tick),  
                      fifelse(!is.na(coop_U) & !is.na(defector_U), (coop_U + defector_U) /2, 
                      fifelse(is.na(coop_U), defector_U, coop_U))]
out[,'av_utility'] <- sum(output[-1, fifelse(!is.na(coop_U) & !is.na(defector_U), (coop_U + defector_U),
                                fifelse(is.na(coop_U), defector_U, coop_U))], 
                      na.rm=TRUE) / output[which.max(tick), tick]


# correct ticked for times when stopping rule was applied! 
if (length(c(                                                # avoid error with numeric(0) 
          tail(output[which(tie_updated == 1),tick], n=1),   # only occurs if neither a tie or strat was updated... 
          tail(output[which(strat_updated == 1),tick], n=1)) # which is apparently a thing
        ) > 0 ) {
          out[, 'ticked' := fifelse(
            output[which.max(tick), tick] !=   
            (max(c(                                      # avoid error with numeric(0) 
                  tail(output[which(tie_updated == 1),tick], n=1),
                  tail(output[which(strat_updated == 1),tick], n=1))  
            )),
            (max(c(                                      # avoid error with numeric(0)  
                  tail(output[which(tie_updated == 1),tick], n=1),
                  tail(output[which(strat_updated == 1),tick], n=1)))
            ),
            ticked)
            ]
        }

if(version2 == "actually_rep") {
    # count strat swap path enters (only relevant for actually_recplicating_min)
    out[, 'num_strat_path_entered'] <- sum(output[strat_path_entered == 1, 
                                        strat_path_entered]) 
    # print(paste0("sanity check, here is the num_strat_path_entered!", out[, 'num_strat_path_entered']))            
}

# process networks
start <- as.matrix(network_list[[1]])         # get initial network 
strats <- node_list[tick == 0, strategy]      # get strats for each agent
prop_ties[1, 'cc'] <- sum(start[strats == 1, strats == 1])/2 # C-C ties 
prop_ties[1, 'cd'] <- sum(start[strats == 1, strats == 0])   # C-D ties
prop_ties[1, 'dd'] <- sum(start[strats == 0, strats == 0])/2 # D-D ties at t0
rm(start, node_list)
for (i in 2:length(network_list)) {   # for all other ticks 
  if (is.null(network_list[[i]]) == FALSE) {    # that have breakdown stored
    prop_ties[i + 1, 'cc'] <- network_list[[i]][[1]]
    prop_ties[i + 1, 'cd'] <- network_list[[i]][[2]]
    prop_ties[i + 1, 'dd'] <- network_list[[i]][[3]]
    }
}
prop_ties[prop_ties == 'NaN' ] <- NA         # fill in gaps
total <- sum(prop_ties[1,])                  # and convert counts to props.
prop_ties[, cc := cc[1], by = cumsum(!is.na(cc))][, cc:= cc / total]
prop_ties[, cd := cd[1], by = cumsum(!is.na(cd))][, cd := cd / total]
prop_ties[, dd := dd[1], by = cumsum(!is.na(dd))][, dd := dd / total]
prop_ties[, time := .I - 1]
# get proportions at final tick 
out[, 'final_CC']  <- prop_ties[time ==  out[,ticked], cc]
out[, 'final_CD']  <- prop_ties[time ==  out[,ticked], cd]
out[, 'final_DD']  <- prop_ties[time ==  out[,ticked], dd]

# # store summary plots of:
# # 1. the proportion of tie types over time 
# prop_ties_time <- ggplot(prop_ties[dd < 1, ], aes(x = time)) +  
#   geom_line(aes(y = cc), color = "green") + 
#   geom_line(aes(y = cd), color = "orange") + 
#   geom_line(aes(y = dd), color = "red") + 
#   theme_classic() + ylab("proportion of ties") + 
#   ggtitle(paste0("rep ", rep, ", degree ", deg,
#                  ", rho ", rho, ", resource inflow ", ri))
# prop_ties_time
# # 2. proportion of agents cooperating 
# prop_c_time <- ggplot(output, aes(x = tick)) +        
#   geom_line(aes(y = prop_coop), color = "green") +     
#   ylab("proportion cooperating") + theme_classic() + 
#   ggtitle(paste0("rep ", rep, ", degree ", deg,
#                  ", rho ", rho, ", resource inflow ", ri))
# prop_c_time

# ------------------------------ EXPORT -------------------------------------#
name <- paste0(paste("smoosh",version2,network_dynamics,ostracism_type,
                     "degree", deg, "ri", ri, "rho", rho,"rep", rep,
                     sep = "_"), ".RData")
save(out, prop_ties, file = name)
     # prop_ties_time, prop_c_time,
    
print(name)
rm(name)



