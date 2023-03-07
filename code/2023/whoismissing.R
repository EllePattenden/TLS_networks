##-------------------------------------------------------------------------#
# WHO DIDN'T FINISH ?
##-------------------------------------------------------------------------#

#Rscript whoismissing.R

suppressMessages(library(data.table))
suppressMessages(library(tidyverse))

# loads as combo
load("/data/gpfs/projects/punim1783/output/floyds_mix/test_combo_chunked.Rda") 

for (a in 1:length(combo)) { 
  combo[[a]][,"index" := .I]    # add index col, will indicate which runs 
}                               # are missing at end. 

finished.list <- list.files("/data/gpfs/projects/punim1783/output/2023/raw_output/rho_sep/NIM/original", 
                        recursive = FALSE,    # no sub folders 
                        pattern = ".Rdata",   
                        full.names = FALSE    # only need strings 
                        )

for (f in finished.list) {
  # Identify the run! 
  # file names: "actually_replicating_min_rep_XX_degree_XX_RI_XX_rho_XX_.Rdata"
  # get the XXs and store 
  RI <- str_match(f, "RI_\\s*(.*?)\\s*_rho")[,2]
  Rho <- str_match(f, "_rho_\\s*(.*?)\\s*_.R")[,2]
  deg <- str_match(f, "degree_\\s*(.*?)\\s*_RI")[,2]
  rep <- str_match(f, "rep_\\s*(.*?)\\s*_degree")[,2]

  print(paste("file #", match(f, finished.list),"degree", deg, "| ri", RI, "| rho", Rho,"| rep", rep, sep = " "))
  
  target <- 69696969
  out <- FALSE 
  attempts <- 1
  
  while (out == FALSE & attempts <= length(combo)) { 
    for (a in 1:length(combo)) { 
      dat <- combo[[a]]
      row <- dat[replicate == rep & degree == deg & resource_inflow == RI & rho == Rho, 
                 index]      # get row # where conditions are met
      if (length(row) != 0) {       # row will = integer(0) if match not found 
        target <- row 
        #dat[replicate == rep & degree == deg & resource_inflow == RI & rho == Rho, index]
        combo[[a]] <- dat[index != target] # drop row using index col
        #print(paste0("chunk=", a, " row=", target))
        out <- TRUE
      }
      #print(paste0("trying chunk", a))
      attempts <- attempts + 1
    }
  }
}
save(combo, file = "/data/gpfs/projects/punim1783/output/2023/raw_output/rho_sep/NIM/original/rerun.Rda")

for (a in 1:length(combo)) { 
  print(gsub(" ", "", toString(dput(as.numeric(combo[[9]]$index))), fixed = TRUE))  # print out so can keep seed (for reproducibility)
}
