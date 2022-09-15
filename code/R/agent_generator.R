
# Create a population of agents             
agent_generator <- function(N,               # population size       
                            initial_coop,    # initial proportion cooperating
                            optimal_effort,  # optimal group effort
                            mu) {            # effort multiplier for defectors 
  agents <- setDT(data.frame(id = 1:N,       # tracking each agent's id,
                             strategy = 0,   # strategy (1 = cooperate)
                             effort = 0,     # effort level
                             payoff = 0,     # payoff 
                             utility = 0     # utility 
  ))
  set(agents,                                 # initialise strategies by  
      sample(nrow(agents), initial_coop * N), # randomly setting cooperators
      "strategy", 1)                          # strategy ==  1
  e_coop <- optimal_effort / N                # set efforts 
  e_defect <- e_coop * mu
  agents[, effort := fifelse(strategy == 0, e_defect, e_coop)]
  return(agents)
}

# Create initial network structure
network_generator <- function(network_type, N, degree) {
  if (network_type == "ER_random") {
    network <- as_adjacency_matrix(
      igraph:::sample_gnm(          # using the "gnm" method
        N,                          # number of nodes
        (N * degree) / 2,           # total number of ties
        directed = FALSE),          # ( 'degree' = average degree for all nodes )
      sparse = TRUE       # get out N x N sparse Matrix of class "dgCMatrix"
    )
  }
  return(network)
}
