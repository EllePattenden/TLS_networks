get_tie_strat <- function(network_dynamics) {
  
  if (network_dynamics == "min_method") {
    return(
      function(agents_own, comp, network) {
          net <- network
          net[comp[agent == "i", id],           # agent_ii severs the tie 
              comp[agent == "j", id]] <- FALSE  # noting this is regardless of utility diff!
          net[comp[agent == "j", id],           # agent_jj also needs to sever the tie... 
              comp[agent == "i", id]] <- FALSE  # because we're working with un-directed ties
          
          # identify  new partner for agent_ii
          options <- setdiff(seq(1, 100, 1),                 # identify potential new partners
                             c(comp[agent == "i", id],       # not self
                               comp[agent == "j", id],       # agent_jj 
                               which(net[comp[agent == "i", id], ] == TRUE) # or current ties
                             ))
          if (length(options) != 0 ) {   # avoiding error in k<90> with SUPER popular agents... 
            new_mate <- sample(options, 1)         
            net[comp[agent == "i", id],       # make new tie
                    new_mate] <- TRUE             # selected randomly from options
            net[new_mate, comp[agent == "i", id]] <- TRUE  
            net[comp[agent == "i", id], new_mate] <- TRUE 
            rm(new_mate)
          }
          return(net)
        }
      )
    }  # end "min_method"
  
  if (network_dynamics == "NIM") {
    return(
      function(agents_own, comp, network) {
          net<- network
          options <-  # find agent_ii's cooperative partners who aren't tied to agent_jj 
            setdiff(agents_own[id %in% which(net[comp[agent == "i", id], ] == TRUE) & strategy == 1, id], 
                    agents_own[id %in% which(net[comp[agent == "j", id], ] == TRUE), id])
          
          if ( # as long as ii has a cooperative partner to recruit 
            length(options) > 0) { 
            agent_kk <- sample(options, 1)  # pick one randomly; they are agent_kk
            
            print(paste0("recruited_mate is", agent_kk)) # sanity check! 
            
            if (       # and as long as kk has the potential to join
              length(  # (i.e., more ties than just to ii)
                which(net[agent_kk, ] == TRUE)[which(net[agent_kk, ] == TRUE) !=  
                                                 comp[agent == "i", id]]) > 0 ) { 
                  k_to_cull <-         # agent_kk then picks one of their social partners, randomly 
                  # we *could* force this to be a cooperator, 
                  # but it seems reasonable to keep the random random random ... 
                    sample(which(net[agent_kk, ] == TRUE)[
                          which(net[agent_kk, ] == TRUE) != comp[agent == "i", id]], 1) 
                # and cuts their tie to them, to maintain <k>, and makes a tie to agent_jj 
                  if (length(k_to_cull) > 0) { 
                    net[agent_kk, k_to_cull] <- FALSE
                    net[k_to_cull, agent_kk] <- FALSE     
                    net[agent_kk, comp[agent == "j",id]] <- TRUE
                    net[comp[agent == "j",id], agent_kk ] <- TRUE
                  }
              }
          }
          return(net)
      }
    )
  }   # end "NIM"
  
}