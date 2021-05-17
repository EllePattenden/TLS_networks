################## one run replicating Min et al. ################## 

one_run <- function(     # reading in row of values in sequence from combo
  replicate, 
  degree, 
  resource_inflow, 
  rho, 
  N, 
  initi_R, 
  max_R, 
  d, 
  E_opt, 
  q, 
  w,
  alpha, 
  beta, 
  gamma, 
  h, 
  t, 
  g, 
  max_t, 
  p_coop, 
  mu, 
  negative_utilities, 
  update_prob, 
  network_type, 
  scheduling_dynamics, 
  tie_strat, 
  e_coop, 
  e_defect
) {   
  
  # get population of agents 
  agents_own <- agent_generator(N, p_coop, E_opt, mu)
  # and their initial social network 
  network <- network_generator(network_type, N, degree)
  
  # working variables
  c <- resource_inflow   # stick to calling this c 
  R <- initi_R    # current size of the resource    
  E_total <- agents_own[, sum(effort)]   # group's total effort 
  p_coop_i <-  agents_own[, sum(strategy)] / N   # proportion agents 
  # cooperating now
  repeats <- update_prob * N    # agents (potentially) comparing strategies in 
  # each tick; defaults to 1
  
  print_ticks <- c(1, seq(0, max_t, by = 1000))  # used for storing data every "by" rounds 

  # functions
  cobb_douglas <- 
    function(E_total, R, alpha, beta, gamma) {           #  Cobb-Douglas
      product <- gamma * (E_total ^ alpha) * (R ^ beta)  #  for production
      return(product)
    }
  resource_dynamics <- 
    function(R, c, d, max_R, q, E_total) {     # Resource dynamics 
      new <- R + (c - ( d * ( (R / max_R ) ^ 2 ) )) - (q * E_total * R ) 
      return(new) 
    }
  ostracism <- 
    function (h, t, g, p_coop_i) {        # Gompertz for ostracism
      o <- h * exp(t * exp(g * p_coop_i))       
      return(o)
    } 
  
  # for storing output 
  network_d <- array(0, dim = c(N, N, length(print_ticks)))
  node_list <- agents_own[, c("id", "strategy")][, tick := 0]
  output <- setDT(data.frame(matrix(0, max_t, 9)))
  setnames(
    output, new = c(   # with cols: 
      "prop_coop",     # (1) proportion of cooperators,
      "R_now",         # (2) resource size after harvesting,
      "magnitude_ostracism",  # (3) ostracism cost
      "coop_U",        # (4) cooperator utility / payoff
      "defector_pi",   # (5) defector payoff
      "defector_U",    # (6) defector utility,
      "total_e",       # (7) total effort,
      "total_harvest",  # (8) total product for round,
      "strategy_updated"  # (9) strategy updated? 1 = true, 0 = false 
    )) 
  
  # Go procedure 
  for (tick in 1:max_t) {
    
    # calculate current proportion of agents cooperating 
    p_coop_i <-  agents_own[, sum(strategy)] / N 
    
    # calculate the total group effort
    E_total <- agents_own[, sum(effort)]
    
    # update the state of the resource to reflect this 
    R <- resource_dynamics(R, c, d, max_R, q, E_total)  
    
    # calculate the groups total product
    total_product <- cobb_douglas(E_total, R, alpha, beta, gamma)
    
    # calculate individual payoffs 
    agents_own[, payoff := 
                 ((effort / E_total) * total_product) - ( w * effort)]
    
    # store utility for cooperators (same as payoff)
    agents_own[strategy == 1, utility := payoff]
    
    # calculate utility for defectors 
    defectors <- agents_own[strategy == 0, id]  # locate defectors
    ostracism_amount <- vector(0, mode = "numeric")
    if (length(defectors) != 0 && 
        length(defectors) != N ) {   # only run if pop isn't at equ.
      payoffs <- unique(agents_own[, payoff, by = strategy])  # extract payoff for each strat.
      defector_adv <- ( payoffs[strategy == 0, payoff] - 
                          payoffs[strategy == 1, payoff] ) / payoffs[strategy == 0, payoff]
      
      ostracism_amount <- numeric(length(defectors))   
      util <- 0
      for (d in defectors) {                       # loop through and 
        partners <- which(network[d, ] == 1)       # find partners 
        n_coop <- nrow(agents_own[id %in% partners & strategy == 1]) # count number cooperating
        ostracism_i <- ostracism(                  # ostracism is det. by the 
          # h, t, g, n_coop                         # num. of partners cooperating
          # h, t, g, (n_coop / N)           # prop. of pop directely tied to who are coop-ing
          h, t, g, (n_coop / length(partners))    # prop of partners cooperating
          # h, t, g, (n_coop / degree)    # prop of partners cooperating
        )
        util <- agents_own[d, payoff] - (ostracism_i * defector_adv)   # track utility
        agents_own[d, utility := util]             # update
        ostracism_amount[d] <- (ostracism_i * defector_adv)    # add to runninng tally
      }
    }
    
    if (negative_utilities == FALSE) {                           # get rid of 
      agents_own[, utility := fifelse(utility < 0, 0, utility)]  #  negatives 
    }
    
    # strategy updating (+ rewiring for scheduling_dynamics == "min_method")
    strat_switched <- vector(length = repeats, mode = "numeric")
    for (r in 1:repeats) {
      agent_i <- sample(nrow(agents_own), 1)     # randomly select focal agent i 
      strat_i <- agents_own[agent_i, strategy]   # record their current strategy
      mates <- which(network[agent_i, ] == 1)  # get list of possible partners 
      if (length(mates) != 0) {    # as long as i >= 1 tie 
        agent_j <- sample(mates, 1)              # pick one 
        strat_j <- agents_own[agent_j, strategy]   # record their strat
        if (strat_i != strat_j) {
          utility_i <- agents_own[agent_i, utility]  # get utilities
          utility_j <- agents_own[agent_j, utility]
          if (utility_i != 0 || utility_j != 0) {    # avoid errors / 0
            utility_diff <-  (utility_i - 
                                utility_j ) / (abs(utility_i) + abs(utility_j))
            
            if (scheduling_dynamics == "min_method") {
              # update strategy with prob rho, otherwise rewire network (p = 1 - rho)
              if (runif(1,0,1) < rho) {      # decide which route to go down  
                if (utility_diff < 0 && runif(1,0,1) < abs(utility_diff) ) {  
                  agents_own[agent_i, strategy := strat_j]             # strategy updating
                  agents_own[strategy == 0, effort := e_defect]
                  agents_own[strategy == 1, effort := e_coop]
                  strat_switched[r] <- 1
                } 
              } else {   # tie rewiring
                if (strat_i == 1 && 
                    strat_j == 0) {   # only happens when i coops and j defects
                  network[agent_i, agent_j] <- 0      # sever tie
                  options <- setdiff(seq(1, 100, 1), 
                                     c(agent_i, agent_j, mates)) # potentials
                  new_tie <- sample(options, 1)           # make tie randomly
                  network[agent_i, new_tie] <- 1
                }
              }
            }
            if (scheduling_dynamics == "seperate") {
              if (utility_diff < 0 && 
                  runif(1,0,1) < abs(utility_diff) ) {
                agents_own[agent_i, strategy := strat_j]        # strategy updating 
                agents_own[strategy == 0, effort := e_defect]
                agents_own[strategy == 1, effort := e_coop]
                strat_switched[r] <- 1
                p_coop_i <- agents_own[, sum(strategy)] / N 
              }
            }
          }
        }
        rm(list = c("agent_i", "strat_i", "mates", "agent_j", "strat_j"))
      }
    }
    
    # rewiring (when scheduling_dynamics == "seperate")
    if (scheduling_dynamics == "seperate") {
      if (runif(1,0,1) < rho) { 
        agent_i <- sample(nrow(agents_own), 1)     
        strat_i <- agents_own[agent_i, strategy]   
        matez <- which(network[agent_i, ] == 1)
        if (length(matez != 0)) {
          agent_j <- sample(matez, 1)                
          strat_j <- agents_own[agent_j, strategy]   
          if (tie_strat == "min_method" )  {
            if (strat_i == 1 && strat_j == 0) {
              network[agent_i, agent_j] <- 0      # sever tie
              options <- setdiff(seq(1, 100, 1), 
                                 c(agent_i, agent_j, matez)) 
              new_tie <- sample(options, 1)       # get new connection
              network[agent_i, new_tie] <- 1
            }
          }
        }
      }
    }
    
    # Stopping rules (and read-out for debugging)
    if (p_coop_i == 1) {
      # print(paste0("Cooperators win at tick", tick, "!"))   # stopping rule 1
      # print("..............................................")
      break
    }
    if (p_coop_i == 0) {
      # print(paste0("defectors win at tick ", tick, "!"))    # stopping rule 2
      # print("..............................................")
      break
    }
    
    
    if (tick %in% print_ticks) {
      nodes <- agents_own[, c("id", "strategy")][, tick := tick]
      l <- list(node_list, nodes)
      node_list <- rbindlist(l, use.names = TRUE)
      network_d[, , match(tick, print_ticks)] <- network
      
      # defectors isolated stopping rule 
      defectors <- agents_own[strategy == 0, id]
      for (d in defectors) {
        partners <- which(network[d, ] == 1)       # find partners 
        n_coop <- nrow(agents_own[id %in% partners & strategy == 1]) # count number cooperating
        if (n_coop == 0) {break} # if none are, stop 
      }
    }
    
    # store main data for round 
    output[tick, names(output) := .(p_coop_i, R, 
                                    mean(ostracism_amount, na.rm = TRUE), 
                                    agents_own[strategy == 1, mean(payoff)],
                                    agents_own[strategy == 0, mean(payoff)],
                                    agents_own[strategy == 0, mean(utility)],
                                    E_total, total_product, sum(strat_switched)) ] 
    rm(list = c("defectors", "strat_switched", "total_product", "E_total"))
  }
  return(list(output, network_d, node_list))
}