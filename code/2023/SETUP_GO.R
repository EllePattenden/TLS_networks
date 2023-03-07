##---------------------------------------------------------------------------
##                          FUNCTIONS FOR GO                               --
##                replicating/refuting/extending Min et al.                --
##---------------------------------------------------------------------------

one_run <- function(        # reading in row of values in order from combo
  replicate, degree, resource_inflow, rho, N, initi_R, max_R, d, E_opt, 
  q, w, alpha, beta, gamma, h, t, g, max_t, p_coop, mu, negative_utilities, 
  update_prob, network_type, e_coop, e_defect, version, network_dynamics, 
  ostracism_type
){
  
  # --------------------------- Set-Up  ------------------------------------#
  
  # Get population of agents 
  agents_own <- agent_generator(N, p_coop, E_opt, mu)
  # and their initial social network (ER-random with given average degree)
  network <- network_generator(network_type, N, degree)
  
  # Initialise working variables 
  c <- resource_inflow      # stick to calling this c for convenience
  R <- initi_R              # current size of the resource
  community_stopping <-     # isolated defectors stopping rule 
    seq(1000, max_t, 1000)  # checked every 1000 ticks 
  # repeats <- update_prob * N  # number selected for each social interaction
  
  # Get sub-procedures (functions): 
  
  ## Cobb-Douglas Production 
  cobb_douglas <- 
    function(E_total, R, alpha, beta, gamma) {          
      product <- gamma * (E_total ^ alpha) * (R ^ beta)  
      return(product)
    }
  
  ## Resource Dynamics 
  resource_dynamics <- 
    function(R, c, d, max_R, q, E_total) {     
      new <- R + (c - ( d * ( (R / max_R ) ^ 2 ) )) - (q * E_total * R ) 
      if (new < 0) { new <- 0 }     # to avoid utility_diff error  
      return(new) 
    }
  
  ## Gompertz for Ostracism 
  ostracism <- 
    function (h, t, g, 
              n_c     # ostracism_type determines how n_c is calculated
              ) {       
      o <- h * exp(t * exp(g * n_c))   
      return(o)
    } 
  
  ## How n_c for ostracism function is calculated 
  calc_nc <- 
    get_nc(ostracism_type)  # sourced from ostracism_functions.R
  
  ## Network Dynamics 
  tie_updating <- 
    get_tie_strat(network_dynamics) # sourced from network_dynamics_functions.R 
  
  # Initialise output storage 
  
  # Networks stored in network_list 
  # Note - we only store network data at t0 and in rounds where ties or strategies
  # are updated. Empty [[]] are filled in later on. 
  network_list <- vector(mode = "list", length = max_t + 1)  
  # There are two options for what data gets chucked in network_list: 
  # 1 = "Big Data", where info about ties are stored as sparse adjacency matrices
  # 2 = "min_method", stores the proportion of CC, CD, and DD ties 
  # 2 is preferable by default, because we can produce a fuck ton of networks in max_t...
  # but you can flip comments in GO if you want 1 
  
  node_list <- data.table(
    id = rep(seq(1, 100, 1), times = 1 + max_t),  # info about key attribute 
    strategy = NA,                                # = strategy
    tick = rep(c(0, 1:max_t), each= N)            # stored in node_list
  ) # will remove rows with strategy == NA at end
  
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
  output[, names(output) := lapply(.SD, as.numeric)]   # make cols. numeric
  
  # store at t0
  output[1 , 
         names(output) := .(        # only relevant values
              (agents_own[, sum(strategy)] / N),   # p_coop_i = .99
              R,                                   # initial size of the resource
              NA, NA, NA, NA, NA, 
              agents_own[, sum(effort)],           # E_total with 99% cooperators 
              NA, NA, NA
              )
    ]
  network_list[[1]] <- network            # initial network
  nodes <- agents_own[, c("id", "strategy")][, tick := 0]
  node_list[nodes, on = c('id', 'tick'), 
            strategy := i.strategy]         # initial strategies 
  
  
  # ------------------------------- GO! ------------------------------------#
  
  for (tick in 1:max_t) { 
    
    # print(paste0("round ", tick))  # for debugging 
    
    # set up 
    strat_switched <- 0 
    tie_switched  <- 0       
    
    # calculate current proportion of agents cooperating 
    p_coop_i <-  agents_own[, sum(strategy)] / N 
    
    # calculate the total group effort
    E_total <- agents_own[, sum(effort)]
    
    # update the state of the resource to reflect this (harvest)
    R <- resource_dynamics(R, c, d, max_R, q, E_total)      
    # 28.09.22 - this should really be at the end of GO ..? 
    # but it's only impacting the first tick so we'll leave it... 
    
    # calculate the groups total product
    total_product <- cobb_douglas(E_total, R, alpha, beta, gamma)
    
    # calculate individual payoffs 
    agents_own[, payoff := 
                 ((effort / E_total) * total_product) - ( w * effort)]
    
    # store utility for cooperators (same as payoff)
    agents_own[strategy == 1, utility := payoff]
    
    # calculate utility for defectors 
    if ( agents_own[, .N, strategy][,.N] == 2) {    # if both strategies are present   
      defector_adv <-  
        agents_own[ , payoff[1], by = strategy ][  # get info for each agent type
          order(-strategy) ][      # make defectors second row  
        ,(V1 - V1[.I - 1]) / V1    # calc. their payoff advantage
      ][2]                     # extract correct (2nd) value
      # locate defectors 
      defectors <- agents_own[strategy == 0, id]
      ostracism_amount <- numeric(length(defectors))  # set up to store
      # loop through each 
      for (d in defectors) {
        util <- agents_own[d, payoff]    # setup util with default value for cases 
                                         # where d doesn't have any ties...  
        if (length(which(network[d, ] == TRUE)) > 0 ){
          n_c <- calc_nc(agents_own, network, d)    # method determined by ostracism_type 
          ostracism_i <- ostracism(h, t, g, n_c)                     
          util <- agents_own[d, payoff] - (ostracism_i * defector_adv)   # track utility
          ostracism_amount[[which(d == defectors)[[1]]]] <- (ostracism_i * defector_adv)   
          rm(ostracism_i, n_c)
        } # else {
          #print(paste0("defector - ", d, " should have no partners below:"))
          #print(partners)
        #}
        agents_own[d, utility := util]  # update
        rm(util)
      }
    } else {  #  if one strategy has dominated,
              # record output and network data then end run
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
    
    # round up negative utilities to 0 (if following TLS, 2016)
    # agents_own[, utility := fifelse(utility < 0, 0, utility)]
    
    # a social interaction, that may culminate in strategy updating, takes place
    agent_i <- agents_own[           # randomly select focal agent i
      sample(nrow(agents_own), 1),   
      .(id, strategy, utility)][     # record their id, strategy, utility
        , agent := "i" ]             # and add an identifier 
                     
    if (length(which(network[agent_i[, id], ] == TRUE)) != 0) { # if agent_i has ties
      agent_j <-                     # randomly select their social partner, agent_j
        agents_own[sample(which(network[agent_i[, id], ] == TRUE), 1),
        .(id, strategy, utility)][
          , agent := "j"]              
      
      comp <- rbind(agent_i, agent_j)  
      rm(agent_i, agent_j)    
      
      if (comp[agent == "i", strategy] !=   # agent_i and agent_j are playing 
          comp[agent == "j", strategy]) {   # different strategies,
        
        # need to avoid error in next step when comp[agent == "i", utility] = 0 
        
        if (comp[agent == "i", utility] == 0) { 
            comp[agent == "i", utility := utility + .0000000001]    # this is another reason
            comp[agent == "j", utility := utility + .0000000001]    # why min et al. annoys me 
          }
        
        utility_diff <-  # agent_i evaluates their utility difference 
          fifelse(                            # calculated as per Min et al. 
            comp[agent == 'i', utility] <     # which is NOT consistent      
            comp[agent == 'j', utility],      # with Shulter et al. (2016)      
            min(1, ((comp[agent == "j", utility] - comp[agent == "i", utility]) / 
                      (abs(comp[agent == "i", utility]) )), 
                 na.rm = TRUE),               # something to come back to...
            0)
    
        # Di in Shulter et al. (2016) 
        #((comp[agent == "j", utility] - comp[agent == "i", utility]) /  
        #          abs(comp[agent == "i", utility]) + abs(comp[agent == "j", utility]))
        # paper has it as Uj = Ui but then that would be a negative number? 
        
        # and updates strat w probability proportional to the payoff difference
        if ( runif(1, min = 0, max = 1) < utility_diff ) {
          agents_own[comp[agent == "i", id], 
                     strategy := comp[agent == "j", strategy]]
          agents_own[comp[agent == "i", id], 
                     effort := fifelse(strategy == 1, e_coop, e_defect)]  
          
          #print(paste("strategy updated!", " Agent", comp[agent == "i",id], 
          # "became a ", agents_own[comp[agent == "i", id], strategy]))
          
          # update trackers 
          strat_switched <- 1
          p_coop_i <- agents_own[, sum(strategy)] / N
          nodes <- agents_own[, c("id", "strategy")][, tick := tick]
          # 23/06/22 changing what network data is stored to limit data size        
          # network_list[[tick + 1]] <- network    # store network
          copy <- as.matrix(network) 
          strat <- nodes[, strategy]
          cc <- sum(copy[strat == 1, strat == 1])/2   # need to divide by 2
          cd <- sum(copy[strat == 1, strat == 0])  
          dd <- sum(copy[strat == 0, strat == 0])/2   # need to divide by 2
          network_list[[tick + 1]] <- list(cc, cd, dd)
          rm(copy, strat, cc, cd, dd)
          # fin changes for storage 
          node_list[nodes, on = c('id', 'tick'), strategy := i.strategy]      # store attributes 
          rm(comp)  # clean up
        }
      }          # end if statement comparing strategies 
    } else {     
      # if agent_i doesn't have ties, they stay that way 
      rm(agent_i)
    }  # end strategy updating procedure 
    
    # social interaction, in which tie updating may happen, transpires
    if (runif(1, 0, 1) < rho) {        # with probability rho 
      agent_ii <- agents_own[
        sample(nrow(agents_own), 1),   # randomly select focal agent ii
        .(id, strategy, utility)][     # record their id, strategy, utility
      , agent := "i" ]                 # and add an identifier
      
      if (  # if agent_ii has ties
        length(which(network[agent_ii[, id], ] == TRUE)) != 0) { 
        agent_jj <- agents_own[          # randomly select one 
          sample(which(network[agent_ii[, id], ] == TRUE), 1), # of their direct neighbors 
          .(id, strategy, utility)       # record " " 
        ][, agent := "j" ]               # and add identifier
        
        comp <- rbind(agent_ii, agent_jj)  
        rm(agent_ii, agent_jj)
        if (comp[agent == "i", strategy == 1] &     # if agent_ii is a cooperator 
            comp[agent == "j", strategy == 0])  {   # and agent_jj is a defector 
                # ties updated, with method determined by "network_dynamics"
                network <- tie_updating(agents_own, comp, network)
                # update trackers 
                tie_switched <- 1 
                nodes <- agents_own[, c("id", "strategy")][, tick := tick]
                #limiting data size by not storing networks, flip comments if desired 
                # network_list[[tick + 1]] <- network    # store network
                copy <- as.matrix(network) 
                strat <- nodes[, strategy]
                cc <- sum(copy[strat == 1, strat == 1])/2   # need to divide by 2
                cd <- sum(copy[strat == 1, strat == 0])  
                dd <- sum(copy[strat == 0, strat == 0])/2   # need to divide by 2
                network_list[[tick + 1]] <- list(cc, cd, dd)
                rm(copy, strat, cc, cd, dd)     # fin changes 
                node_list[nodes, on = c('id', 'tick'), strategy := i.strategy]    # store attributes 
              }
            rm(comp)
            }  # if agent_ii doesn't have ties 
            rm(agent_ii)
          }  # end tie updating procedure 
    
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
        break      # stop 
        }                     
    }
    
    rm(list = c("E_total", "total_product", "strat_switched", "tie_switched", 
                "defectors", "p_coop_i", "nodes"))
    
    if (tick == max_t) {
      print("... Made it to the end! ... ")
    }
    
  }  # end tick 
  
  # --------------------------- Write Data ---------------------------------#
  
  node_list[, strategy := as.numeric(strategy)]  # don't want logical 
  
  name <- paste("NIM_replicating_min", "rep", replicate, 
                "degree", degree, 
                "RI", resource_inflow,
                "rho", rho,
                ".Rdata", sep = "_")
  data <- list(output, network_list, node_list)
  save(data, file = name)
  
  # clean up 
  rm(output, network_list, node_list, # where, this, 
     name, data, c, R)
  
}  # end run 