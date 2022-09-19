plot_time <- function(path, 
                      ri, 
                      rho, 
                      degree,
                      which_plot   # options: "prop_c" & "prop_ties"
                      ) { 
  require("tidyverse")
  require("ggpubr")
  
  files <- list.files(
    as.character(path),
    recursive = FALSE,    
    pattern = paste0("^summary_degree_", degree, "_ri_", ri, "_rho_",
                     rho,"_"), 
    full.names = TRUE
  )
  static <- vector('list',           # name comes from when I 
                   length(files))    # exploring static networks 
  for (f in 1:length(files)) {       # too lazy to change now... 
    load(files[[f]])    # loads 'out', prop_ties, prop_c_time, prop_ties_time
    if (which_plot == "prop_c") {
      static[[f]] <- prop_c_time
    } else {
      static[[f]] <- prop_ties_time
    }
    rm(out, prop_ties, prop_c_time, prop_ties_time)
  }
  hold <- ggarrange(plotlist = static, ncol = 4, nrow = 5)
  assign(paste0("degree_",degree,"_ri",ri,"_rho", rho) , hold)
  rm(hold)
  ggexport(get(paste0("degree_",degree,"_ri",ri,"_rho", rho)), 
           ncol = 4,
           width = 1150, height = 1600,
           pointsize = 12,
           filename = paste0(which_plot,"_degree_",degree,"_ri",ri,"_rho", rho, ".png")
           )
  # return(get(paste0("degree_",degree,"_ri",ri,"_rho", rho)))
  }