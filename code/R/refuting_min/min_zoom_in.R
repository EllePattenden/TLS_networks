##--------------------------------------------------------------------------
##                      REPLICATE/REFUTE MIN ET AL.                       --
##                               zooming in                               --
##--------------------------------------------------------------------------
#install.packages("viridis")  # Install
library(data.table)
library(viridis)           
load("./summary_min.RData")   # loads as summary
summary <- summary[-1, ]

# calc. proportion of ticks where ties and strats were updated 
summary[, prop_ties_switched := num_tie_swaps / ticked]   # ties
summary[, prop_strats_switched := num_strat_swaps / ticked] # strats

# heat maps 
means_wanted <- c(names(summary)[5], names(summary)[14:15])
means <- summary[, lapply(.SD, mean, na.rm =T), 
                 by = c('ri', 'rho', 'degree'), .SDcols = (means_wanted)]
means[, fpc_rounded := round(f_prop_coop, 1)] 

# average proportion cooperating at fin
means %>% filter(degree == 15) %>%
  #ggplot(aes(x = ri, y = rho, fill = f_prop_coop)) +
  ggplot(aes(x = ri, y = rho, fill = fpc_rounded)) +   # rounded like in Min 
  geom_tile() + 
  theme_classic() + 
  scale_x_discrete(breaks = seq(0,60, 10)) +
  scale_y_discrete(breaks = seq(0,1,.1)) +
  labs(fill = "cooperating") +
  xlab("resource inflow") + 
  ylab("rewiring probability") + 
  scale_fill_viridis(option = "plasma") 

means %>% filter(degree == 90) %>%
  #ggplot(aes(x = ri, y = rho, fill = f_prop_coop)) +
  ggplot(aes(x = ri, y = rho, fill = fpc_rounded)) +   # rounded like in Min 
  geom_tile() + 
  theme_classic() + 
  scale_x_discrete(breaks = seq(0,60, 10)) +
  scale_y_discrete(breaks = seq(0,1,.1)) +
  labs(fill = "cooperating") +
  xlab("resource inflow") + 
  ylab("rewiring probability") + 
  scale_fill_viridis(option = "plasma") 



# average prob of strategy updating across run 
# is APPROXIMATE - see stopping rule for isolated defectors 
means %>% filter(degree == 15) %>%
  arrange(ri) %>%
  ggplot(aes(x = ri, y = rho, fill = prop_strats_switched)) +
  geom_tile() + 
  theme_classic() + 
  scale_x_discrete(breaks = seq(0,60, 10)) +
  scale_y_discrete(breaks = seq(0,1,.1)) +
  labs(fill = "strategy updating") +
  xlab("resource inflow") + 
  ylab("rewiring probability") + 
  scale_fill_viridis(option = "mako") 

means %>% filter(degree == 90) %>%
  arrange(ri) %>%
  ggplot(aes(x = ri, y = rho, fill = prop_strats_switched)) +
  geom_tile() + 
  theme_classic() + 
  scale_x_discrete(breaks = seq(0,60, 10)) +
  scale_y_discrete(breaks = seq(0,1,.1)) +
  labs(fill = "strategy updating") +
  xlab("resource inflow") + 
  ylab("rewiring probability") + 
  scale_fill_viridis(option = "mako") 



# ------------------------------------------------------ COME BACK TO 

# geom_area/ribbon
# find xy cooridinates for each prob band 
# method: for each rho, what is the max ri before a shift occurs? 
working <- matrix(0, length(seq(0, 1, .05)), length(seq(0, 1, .1)))
colnames(working) <- c(as.character(seq(0, 1, .1)))   # cols = prob_bands
rownames(working) <- c(as.character(seq(0, 1, .05)))  # rows = rho 
for (i in seq(0,1,.05)) {   # rho starting from 0
  w2<- means[rho == i, c("ri","fpc_rounded")]
  for (j in seq(0, 1,.1)) {  # prob_band 
    working[as.character(i),as.character(j)] <- max(w2[fpc_rounded == j, ri])
  }
}
working <- rownames_to_column(as.data.frame(working), "rho")
# convert to long 
working <- gather(as.data.frame(working), key = "fpc_rounded", value = "max_ri", 2:12)
# giving up on this for now... 

