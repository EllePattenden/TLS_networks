##---------------------------------------------------------------------------
##                          FUNCTIONS FOR GO                               --
##                    replicating/refuting Min et al.                      --
##---------------------------------------------------------------------------

one_run <- function(        # reading in row of values in sequence from combo
  replicate, degree, resource_inflow, rho, N, initi_R, max_R, d, E_opt, q, w,
  alpha, beta, gamma, h, t, g, max_t, p_coop, mu, negative_utilities, 
  update_prob, network_type, scheduling_dynamics, tie_strat, e_coop, e_defect
){
  
  # --------------------------- Set-Up  ------------------------------------#
  # get population of agents 
  agents_own <- agent_generator(N, p_coop, E_opt, mu)
  # and their initial social network (ER-random with given average degree)
  network <- network_generator(network_type, N, degree)
  
  # initialise working variables 
  c <- resource_inflow      # stick to calling this c for convenience
  R <- initi_R              # current size of the resource
  community_stopping <- seq(1000, max_t, 1000)  # isolated defectors stopping  
  # rule checked every 1000 ticks 
  # repeats <- update_prob * N  # number selected for each social interaction
  
  # functions: 
  # Cobb-Douglas production 
  cobb_douglas <- 
    function(E_total, R, alpha, beta, gamma) {          
      product <- gamma * (E_total ^ alpha) * (R ^ beta)  
      return(product)
    }
  # resource dynamics 
  resource_dynamics <- 
    function(R, c, d, max_R, q, E_total) {     
      new <- R + (c - ( d * ( (R / max_R ) ^ 2 ) )) - (q * E_total * R ) 
      if (new < 0) {
        new <- 0        # 23.06.22 - avoid utility_diff error ! 
      }
      return(new) 
    }
  # Gompertz for ostracism 
  ostracism <- 
    function (h, t, g, p_coop) {       
      o <- h * exp(t * exp(g * p_coop))   
      return(o)
    } 
  
  # initialise output storage 
  network_list <- vector(     # note - will only store network data at t0 and in 
    mode = "list",            # rounds where ties or strategies are updated 
    length = max_t + 1        # info about ties stored as sparse adjacency matrices
  )                           # in network_list 
  # (will remove empty [[]] at end)
  node_list <- data.table(
    id = rep(seq(1, 100, 1), times = 1 + max_t),  # info about key attribute 
    strategy = NA,                                # = strategy
    tick = rep(c(0, 1:max_t), each= N)            # stored in node_list
  )                                   # will remove rows with strategy == NA at end
  
  output <- setNames(     # main output from every tick 
    data.table(           # stored as data.table: 
      matrix(nrow = max_t + 1, # output                           
             ncol = 11)),  # with 11 cols: 
    c("prop_coop",       # (1) proportion of cooperators
      "R_now",           # (2) resource size after harvesting,
      "mag_ostracism",   # (3) ostracism cost (average)
      "sd_ostracism",    # (4) ostracism cost (standard deviation)
      "coop_U",          # (5) cooperator utility / payoff  (average*)
      "defector_pi",     # (6) defector payoff (average*)
      "defector_U",      # (7) defector utility (average)
      "total_e",         # (8) total effort,
      "total_harvest",   # (9) total product for round,
      "strat_updated",   # (10) strategy updated? 1 = true, 0 = false 
      "tie_updated"     # (11) tie updated? 1 = true, 0 = false
    ))
  output[, names(output) := lapply(.SD, as.numeric)]   # make all cols. numeric
  
  # store at  t0
  output[1 , names(output) := . (        # only relevant values
    (agents_own[, sum(strategy)] / N),   # p_coop_i = .99
    R,                                   # initial size of the resource
    NA, NA, NA, NA, NA, 
    agents_own[, sum(effort)],           # E_total with 99% cooperators 
    NA, NA, NA)]
  network_list[[1]] <- network              # initial network
  nodes <- agents_own[, c("id", "strategy")][, tick := 0]
  node_list[nodes, on = c('id', 'tick'), 
            strategy := i.strategy]         # initial strategies 
  
  
  # ------------------------------- GO! ------------------------------------#
  
  for (tick in 1:max_t) { 
    
    # set up 
    strat_switched <- 0 
    tie_switched  <- 0       
    
    # calculate current proportion of agents cooperating 
    p_coop_i <-  agents_own[, sum(strategy)] / N 
    
    # calculate the total group effort
    E_total <- agents_own[, sum(effort)]
    
    # update the state of the resource to reflect this (harvest)
    R <- resource_dynamics(R, c, d, max_R, q, E_total)  
    
    # calculate the groups total product
    total_product <- cobb_douglas(E_total, R, alpha, beta, gamma)
    
    # calculate individual payoffs 
    agents_own[, payoff := 
                 ((effort / E_total) * total_product) - ( w * effort)]
    
    # store utility for cooperators (same as payoff)
    agents_own[strategy == 1, utility := payoff]
    
    # calculate utility for defectors 
    if (agents_own[, .N, strategy][,.N] == 2) {      # only run if both strategies are present 
      
      #  && agents_own[strategy == 0, payoff[1]] != 0) {   # and the payoff to defectors isn't 0 
      # this occasionally throws an error (assuming when there are no defectors)
      # even though using && should return false if there is only one strategy in the pop... 
      
      defector_adv <- agents_own[, payoff[1], by = strategy   # get info for each agent type
      ][order(-strategy)           # make defectors second row   
      ][,(V1 - V1[.I - 1]) / V1  # calc. their payoff advantage
      ][2]                     # extract correct (2nd) value
      # locate defectors 
      defectors <- agents_own[strategy == 0, id]
      ostracism_amount <- numeric(length(defectors))  # set up to store
      # loop through and find partners 
      for (d in defectors) {
        partners <- which(network[d, ] == TRUE)    # find partners 
        n_coop <- nrow(agents_own[id %in% partners & 
                                    strategy == 1,]) # count number cooperating
        ostracism_i <- ostracism(                  # ostracism is then determined by the 
          h, t, g, (n_coop / length(partners))     # proportion of partners cooperating
        )
        util <- agents_own[d, payoff] - (ostracism_i * defector_adv)   # track utility
        agents_own[d, utility := util]                                 # update
        ostracism_amount[[which(d == defectors)[[1]]]] <- 
          (ostracism_i * defector_adv)                           # add to running total
        rm(util, partners, n_coop, ostracism_i)
      }
    } else { 
      # if one strategy has dominated, record output and networks then end run
      agents_own[strategy == 0, utility := payoff]
      strat_switched <- NA
      tie_switched <- NA 
      ostracism_amount <- 0 
      output[tick + 1 , names(output) :=    # + 1 so all data is recorded 
               .(p_coop_i,   # proportion of cooperators    
                 R,          # size of resource after harvesting
                 mean(ostracism_amount,      # average ostracism 
                      na.rm = TRUE),         # felt by defectors
                 sd(ostracism_amount,        # sd of ostracsim
                    na.rm = TRUE),           # felt by defectors 
                 agents_own[strategy == 1,   # cooperator 
                            mean(payoff)],   # average payoff 
                 agents_own[strategy == 0,   # defector 
                            mean(payoff)],   # average payoff 
                 agents_own[strategy == 0,   # defector 
                            mean(utility)],  # average utility 
                 E_total,           # total effort
                 total_product,     # total product 
                 strat_switched,    # did someone swap strategies?
                 tie_switched       # did someone swap ties? 
               ) ] 
      network_list[[ tick + 1]] <- network  
      nodes <- agents_own[, c("id", "strategy")][, tick := tick]
      node_list[nodes, on = c('id', 'tick'), strategy := i.strategy]
      break 
    }
    # if you want to round up negative utilities to 0 (following Shultz et al., 2016)
    # agents_own[, utility := fifelse(utility < 0, 0, utility)]
    
    # a social interaction that may culminate in strategy updating takes place
    agent_i <- agents_own[
      sample(nrow(agents_own), 1),   # randomly select focal agent i
      .(id, strategy, utility)       # record their id, strategy, utility
    ][, agent := "i"               # and add an identifier
    ]               
    if ( length(which(network[agent_i[, id], ] == TRUE)) != 0) { # if agent_i has ties
      agent_j <- agents_own[         # randomly select their social partner, agent_j
        sample(which(network[agent_i[, id], ] == TRUE), 1),
        .(id, strategy, utility)     # record " " 
      ][, agent := "j"               # and add identifier
      ]
      comp <- rbind(agent_i, agent_j)  
      rm(agent_i, agent_j)
      # agent i evaluates their utility
      if (comp[agent == "i", strategy] !=    # but only if strats differ
          comp[agent == "j", strategy]) {  
        utility_diff <- fifelse(              # calculated as per Min et al. 
          comp[agent == 'i', utility] <       # which is NOT consistent 
            comp[agent == 'j', utility],      # with Shulter et al. (2016)
          min(1, ((comp[agent == "j", utility] -      # something to come back to...
                     comp[agent == "i", utility]) /  
                    abs(comp[agent == "i", utility])), na.rm = TRUE),
          0)
        
        # and updates strat w probability proportional to the payoff difference
        if (runif(1, min = 0, max = 1) < utility_diff) {
          agents_own[comp[agent == "i", id], 
                     strategy := comp[agent == "j", strategy]]
          agents_own[comp[agent == "i", id],                  # do sequentially, as run 
                     effort := fifelse(strategy == 1, e_coop, e_defect)]  # in parallel
          
          print(paste("strategy updated!", " Agent", comp[agent == "i",id], 
                      "became a ", agents_own[comp[agent == "i", id], strategy]))
          
          # update trackers 
          strat_switched <- 1
          p_coop_i <- agents_own[, sum(strategy)] / N
          network_list[[tick + 1]] <- network    # store network
          nodes <- agents_own[, c("id", "strategy")][, tick := tick]
          node_list[nodes, on = c('id', 'tick'),    
                    strategy := i.strategy]      # store attributes 
          rm(comp, nodes)  # clean up
        }
      }
    } else {
      rm(agent_i)
    }  # end strategy updating procedure 
    
    # social interaction in which tie updating may happen transpires
    if (runif(1, 0, 1) < rho) {        # with probability rho 
      agent_ii <- agents_own[
        sample(nrow(agents_own), 1),   # randomly select focal agent ii
        .(id, strategy, utility)       # record their id, strategy, utility
      ][, agent := "i"               # and add an identifier
      ]
      if (length(which(network[agent_ii[, id], ] == TRUE)) != 0) { # if agent_ii has ties
        agent_jj <- agents_own[          # randomly select one 
          sample(which(network[agent_ii[, id], ] == TRUE), 1), # of their direct neighbors 
          .(id, strategy, utility)       # record " " 
        ][, agent := "j"               # and add identifier
        ]
        comp <- rbind(agent_ii, agent_jj)  
        rm(agent_ii, agent_jj)
        if (comp[agent == "i", strategy == 1] &     # if agent_ii is a cooperator 
            comp[agent == "j", strategy == 0])  {   # and agent_jj is a defector
          network[comp[agent == "i", id],           # agent_ii severs the tie 
                  comp[agent == "j", id]] <- FALSE  # noting this is regardless of utility diff!
          
          ## 22/06/22 - past Elle, you dumbO these are NOT directed ties...
          ## so need to update for both rows!
          
          network[comp[agent == "j", id],           # agent_jj also needs to sever the tie... 
                  comp[agent == "i", id]] <- FALSE
          
          options <- setdiff(seq(1, 100, 1),             # identify potential new partners
                             c(comp[agent == "i", id],   # not self
                               comp[agent == "j", id],       # agent_jj 
                               which(network[comp[agent == "i", id], ] == TRUE) # or current ties
                             ))
          
          ## 22/06/22 - changes so that row j is updated as well 
          new_mate <- sample(options, 1)                ## 22/06/22 - done first 
          network[comp[agent == "i", id],       # make new tie
                  new_mate] <- TRUE             # selected randomly from options
          network[new_mate, comp[agent == "i", id]] <- TRUE  
          rm(new_mate)
          
          # update trackers
          tie_switched <- 1 
          rm(copy, strat, cc, cd, dd)
          nodes <- agents_own[, c("id", "strategy")][, tick := tick]
          
          #23/06/22 changing what network data is stored to limit data size
          # network_list[[tick + 1]] <- network    # store network
          copy <- as.matrix(network) 
          strat <- nodes[, strategy]
          cc <- sum(copy[strat == 1, strat == 1])/2   # need to divide by 2
          cd <- sum(copy[strat == 1, strat == 0])  
          dd <- sum(copy[strat == 0, strat == 0])/2   # need to divide by 2
          network_list[[tick + 1]] <- list(cc, cd, dd)
          
          node_list[nodes, on = c('id', 'tick'),    
                    strategy := i.strategy]      # store attributes 
          rm(options, comp, nodes)  # clean up
        } else {
          rm(comp)
        }
      } else {
        rm(agent_ii)
      }
    }   # end tie updating procedure 
    
    # save data for tick 
    output[tick + 1 , names(output) :=    # + 1 so all data is recorded 
             .(p_coop_i,   # proportion of cooperators    
               R,          # size of resource after harvesting
               mean(ostracism_amount,      # average ostracism 
                    na.rm = TRUE),         # felt by defectors
               sd(ostracism_amount,        # sd of ostracsim
                  na.rm = TRUE),           # felt by defectors
               agents_own[strategy == 1,   # cooperator 
                          mean(payoff)],   # average payoff 
               agents_own[strategy == 0,   # defector 
                          mean(payoff)],   # average payoff 
               agents_own[strategy == 0,   # defector 
                          mean(utility)],  # average utility 
               E_total,           # total effort
               total_product,     # total product 
               strat_switched,    # did someone swap strategies?
               tie_switched       # did someone swap ties?
             ) ] 
    
    # Stopping rules (and read-out for debugging)
    if (p_coop_i == 1) {
      print(paste0("Cooperators win at tick", tick, "!"))   # stopping rule 1
      print("..............................................")
      break
    }
    if (p_coop_i == 0) {
      print(paste0("defectors win at tick ", tick, "!"))    # stopping rule 2
      print("..............................................")
      break
    }
    if (tick %in% community_stopping) {
      defectors <- agents_own[strategy == 0, id]
      n_coop_d <- numeric(length(defectors))
      for (d in defectors) {
        partners <- which(network[d, ] == TRUE) # find partners
        n_coop_d[ which(d == defectors)[[1]]] <- 
          nrow(agents_own[id %in% partners & strategy == 1, ]) # store count cooperating
        rm(partners)
      }
      if (sum(n_coop_d) == 0) {      # if all defectors are in isolated community
        print(paste0("defectors in isolated community at", tick, "!"))  
        print("..............................................")
        break }                      # stop 
    }
    
    rm(list = c("E_total", "total_product", "strat_switched", "tie_switched", 
                "defectors", "p_coop_i"))
    
    
    if (tick == max_t) {
      print("... Made it to the end! ... ")
    }
    
  }  # end tick 
  
  # --------------------------- Write Data ---------------------------------#
  
  node_list[, strategy := as.numeric(strategy)]  # don't want logical 
  
  #where <- getwd()
  #this <- paste(# Sys.Date(), 
  # "rep", replicate, 
  # "degree", degree, 
  # "RI", resource_inflow,
  # "rho", rho, 
  # sep = "_")
  
  # method for local 
  #dir.create(file.path(
  # where, 
  # "replicating_min", this), recursive = TRUE)
  
  # name <- file.path(where, 
  # "replicating_min",  # this, 
  # paste("rep", replicate, 
  #     "degree", degree, 
  #    "RI", resource_inflow,
  #   "rho", rho,
  #   ".Rdata", sep = "_"))
  
  name <- paste("replicating_min", "rep", replicate, 
                "degree", degree, 
                "RI", resource_inflow,
                "rho", rho,
                ".Rdata", sep = "_")
  data <- list(output, network_list, node_list)
  save(data, file = name)
  
  # clean up 
  rm(output, network_list, node_list,
     # where, this, 
     name, data, c, R)
  
}  # end run 

##---------------------------------------------------------------------------
##                            PARAMETER SWEEEEEP                           --
##---------------------------------------------------------------------------

# Function to run ^ for all parameter combinations 
# (will move this to mclappy for parallel)
modelRun <- function(combo) {
  mapply(one_run, combo[,replicate], combo[, degree], 
         combo[, resource_inflow], combo[, rho], combo[, N], 
         combo[, initi_R], combo[, max_R], combo[, d], 
         combo[, E_opt], combo[, q], combo[, w], combo[, alpha],
         combo[, beta], combo[,gamma], combo[, h], combo[, t],
         combo[, g], combo[, max_t], combo[, p_coop], combo[, mu],
         combo[, negative_utilities], combo[,update_prob], 
         combo[, network_type], combo[, scheduling_dynamics], 
         combo[, tie_strat], combo[, e_coop], combo[, e_defect]
  )
}
