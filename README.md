# TLS_networks

This repository contains a series of simulations that add social networks and network dynamics to Schlüter, M., Tavoni, A., & Levin, S. (2016). Robustness of norm-driven cooperation in the commons. Proceedings of the Royal Society B: Biological Sciences, 283(1822), 20152431 (https://royalsocietypublishing.org/doi/full/10.1098/rspb.2015.2431), an ABM of Tavoni, A., Schlüter, M., & Levin, S. (2012). The survival of the conformist: social pressure and renewable resource management. Journal of theoretical biology, 299, 152-161. (https://www.sciencedirect.com/science/article/abs/pii/S002251931100347X)

This builds upon work looking at the impact of static networks in this model, including: 
- Chung, N. N., Chew, L. Y., & Lai, C. H. (2013). Influence of network structure on cooperative dynamics in coupled socio-ecological systems. EPL (Europhysics Letters), 104(2), 28003.
- Sugiarto, H. S., Chung, N. N., Lai, C. H., & Chew, L. Y. (2015). Socioecological regime shifts in the setting of complex social interactions. Physical Review E, 91(6), 062804.
- Sugiarto, H. S., Lansing, J. S., Chung, N. N., Lai, C. H., Cheong, S. A., & Chew, L. Y. (2017). Social cooperation and disharmony in communities mediated through common pool resource exploitation. Physical Review Letters, 118(20), 208301.

The starting point is the lone study (to my knowledge) that incorporated network dynamics: Min, Y., Du, Y., & Jin, C. (2018). The effect of link rewiring on a coevolutionary common pool resource game. Physica A: Statistical Mechanics and its Applications, 512, 935-944. 
- Problematic for reasons outlined in main_min.R; first job is clarifying results when the probability of tie rewiring is seperated from the probability of strategy updating. 
- main_min.R recreates figure three from the paper, with rho impacting tie but not strategy updating. 
- need to follow up with same parameter space with rho impacting probability of both, then prob of tie rewiring held constant and variable strategy updating (level above prob based on utility comparisons). 

After this, will (1) introduce different tie updating strategy that aligns with the risk hypothesis and (2) give ties a new purpose (monitoring + sanctioning) 




