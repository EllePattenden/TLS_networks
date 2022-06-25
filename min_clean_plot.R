##--------------------------------------------------------------------------
##                      REPLICATE/REFUTE MIN ET AL.                       --
##                      Data Cleaning and Analysis                        --
##--------------------------------------------------------------------------

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
igraph_options(sparsematrices = TRUE)  # work with sparse matrices, please

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
  new_name <- stringr::str_sub(old[[1]], start = 117) %>% 
    str_remove(".Rdata") %>% str_remove("/")
  assign(sub(old[[1]], new_name, i), get(i))
  rm(i, old, new_name)
}
# this is not reproducible, but will work for now
rm(list = ls(pattern = '/Users/ellepattenden/Dropbox/PhD/2021/simulations/'))

# ------------------------------ Prelim ------------------------------------#

# set up folder for storing output 
if (!dir.exists("results")) {     
  dir.create("results")         
} 

prop_coop_at_fin <- data.table(ri = rep(seq(10, 60, 10), times = 11), 
                               rho = rep(seq(0, 1, .1), each = 6),
                               prop_coop = NA,
                               ticked = NA)
prop_coop_at_fin[, names(prop_coop_at_fin) := lapply(.SD, as.numeric)]

stuff <- ls()[! ls() %in% "prop_coop_at_fin"]

for (i in stuff) { 
  
  # for saving plots / data
  where <- getwd()
  dir.create(file.path(where, "results", i),  
             recursive = TRUE)
  path <- file.path(where, "results", i)
  
  # identify condition
  ri1 <- str_match(i, "RI_\\s*(.*?)\\s*_rho")
  ri1 <- ri1[, 2]
  rho1 <- str_match(i, "_rho_\\s*(.*?)\\s*_")
  rho1 <- rho1[, 2]
  # get main output
  dat <- get(i)[[1]]
  dat[, tick := .I - 1]
  dat <- na.omit(dat, "prop_coop")
  # get ties
  ties <- get(i)[[2]]
  names(ties) <- seq_along(ties) - 1 # label with tick
  ties[sapply(ties, is.null)] <- NULL   # remove empty
  # get attributes
  att <- get(i)[[3]]
  # for plotting
  networks_to_plot <- names(ties)  
  networks_to_plot <- c(        # only eyeball initial 
    head(networks_to_plot, 1),  # and final networks 
    tail(networks_to_plot, 1))  # for now 
  
  # extract prop_coop and tick at fin 
  if (dat[which.max(prop_coop), prop_coop] != 1) {
    # min_coop <- dat[which.min(prop_coop), prop_coop]    # 25/06/22 - this isn't nec. the final coop
    min_coop <- dat[which.max(tick), prop_coop]           # this is 
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
  
  # prelim plots 
  prop_coop_time <- ggplot(dat, aes(x = tick)) +  
    geom_line(aes(y = prop_coop)) +     # proportion cooperating 
    labs(title = as.character(i)) +     # over ticks 
    ylab("proportion cooperating")
  ggsave(filename = "prop_coop_time.png", 
         plot = prop_coop_time, path =  path)
  
  for (net in networks_to_plot) {
    edge_list <- att %>% filter(tick == net) %>% select(id, strategy)
    network <- graph.adjacency(ties[as.character(net)][[1]],
                               mode = "undirected")
    V(network)$strategy <-  edge_list$strategy
    V(network)$color <- ifelse(V(network)$strategy == 1, "green", "red")
    set.seed(1234)
    # plot and save 
    png <- file.path(path, paste0(net, ".png"))
    png(png)
    plot(network, 
         layout = layout_with_lgl, 
         vertex.size = 10, 
         vertex.label.cex = 0.5)
    dev.off()
    rm(png)
  }
  
  # save clean 
  name <- file.path(where, "results", i,paste( i, ".Rdata", sep = ""))
  save(dat, ties, att, file = name)
  
  rm(dat, att, ties, ri1, rho1, where, name, path)
}
rm(stuff, i, edge_list, network, networks_to_plot, net)

# only working with limited number of runs at the moment
prop_coop_at_fin <- na.omit(prop_coop_at_fin, "ticked")

# make heatmap 
figure_three <- 
  ggplot(prop_coop_at_fin, aes(x = ri, y = rho, fill = prop_coop)) + 
  geom_tile() + 
  xlab("resource inflow, c") + 
  ylab("rewiring probability, rho")
ggsave("figure_three.png", path = file.path(getwd(), "results"))


