library(data.table)

args <- commandArgs(trailingOnly=TRUE) 
print(args) 
version <- args[1]
go <- args[2]

files <- list.files(
    if(version == "actually_rep") {
      if (go == "original") {
         "/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/actually_replicating_min/to_smoosh"
      } else {
        "/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/actually_replicating_min/to_smoosh_R50/"
      }   
    } else {
        if (go == "original") {
          "/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/to_smoosh2"    # to_smoosh2 for rho_sep 
      } else {
        "/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/to_smoosh_R50/"
      }   
    }, 
    recursive = FALSE,      # no sub folders 
    pattern = "^smoosh",   
    full.names = TRUE       # need the file path to load
)
#print(files)

summary <- data.table(rep = NA,               # output for ID + plotting 
                      ri = NA,                # resource inflow
                      rho = NA,               # prob. strategy updating
                      degree = NA,            # average k: 15 or 90 
                      f_prop_coop = NA,       # final proportion of cooperators
                      ticked = NA,            # when the run ended
                      av_prop_coop = NA,      # mean prop coop. over the run
                      sd_prop_coop = NA,      # sd of ^ 
                      num_tie_swaps = NA,     # tally times ties swapped
                      num_strat_swaps = NA,    # tally times strategies swapped,
                      final_CC = NA,          # proportion C-C ties at t_fin
                      final_CD = NA,          # proportion C-D ties at t_fin
                      final_DD = NA,           # proportion D-D ties at t_fin
                      num_strat_path_entered = NA
)

for (f in files) { 
  load(f)    # loads as 'out'
  print(f)
  print(out)
  summary <- rbind(summary, out, fill=TRUE)
  rm(out)
}
rm(f)


if(version == "actually_rep") {
    if (go == "original") {
        save(summary, file = "summary_min.RData")
    } else {
        save(summary, file = "summary_min_R50.RData")
    }   
} else {
    if (go == "original") {
        save(summary, file = "RHO_summary_min.RData") 
    } else {
        save(summary, file = "RHO_summary_min_R50.RData") 
    }
}   
    
# make v rough heatmaps
# library(tidyverse)
# summary <- as.data.frame(summary) 
# fig3_15 <- summary %>% 
#   filter(degree == 15) %>% 
#   ggplot(aes(x = ri, y = rho, fill = f_prop_coop)) +
#   geom_tile() + 
#   theme_classic() + 
#   xlab("resource inflow, c") + 
#   ylab("rewiring probability")
# ggsave("figure_three_15.png")
# fig3_90 <- summary %>% 
#   filter(degree == 90) %>% 
#   ggplot(aes(x = ri, y = rho, fill = f_prop_coop)) +
#   geom_tile() + 
#   theme_classic() + 
#   xlab("resource inflow, c") + 
#   ylab("rewiring probability") 
# ggsave("figure_three_90.png")
