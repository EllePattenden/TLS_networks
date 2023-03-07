get_nc <- function(ostracism_type) {
  if (ostracism_type == "original") {
    return(
      function(agents_own, network, agent_d) {
        n_all <- sum(as.matrix(network[agent_d,]))
        n_coop <- agents_own[,strategy] %*% as.matrix(network[agent_d,]) 
        n_c <- n_coop / n_all 
        return(n_c)
        }
    )
    }
  if (ostracism_type == "FEJ") {
    return(
      function(agents_own, network, agent_d) {
        # working on this! 
      }
    )
  }
}