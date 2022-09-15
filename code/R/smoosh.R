library(data.table)

files <- list.files(
  "/data/gpfs/projects/punim1783/output/floyds_mix/cleaned_sim_results/to_smoosh",
  recursive = FALSE,     # don't go looking in results sub 
  pattern = ".RData", 
  full.names = TRUE
  )

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
                      final_DD = NA           # proportion D-D ties at t_fin
)

for (f in files) { 
  load(f)    # loads as 'out'
  print(out)
  summary <- rbind(summary, out)
  rm(out)
  }

save(summary, file = "summary_min.RData")

# make heatmap
library(tidyverse)
summary <- as.data.frame(summary) 
fig3_15 <- summary %>% 
  filter(degree = 15) %>% 
  ggplot(aes(x = ri, y = rho, fill = prop_coop)) +
  theme_classic() + 
  xlab("resource inflow, c") + 
  ylab("rewiring probability")
ggsave("figure_three_15.png")
fig3_90 <- summary %>% 
  filter(degree = 90) %>% 
  ggplot(aes(x = ri, y = rho, fill = prop_coop)) +
  theme_classic() + 
  xlab("resource inflow, c") + 
  ylab("rewiring probability") 
ggsave("figure_three_90.png")


