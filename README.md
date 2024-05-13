# TLS_networks

This repository contains a series of simulations that add social networks and network dynamics to [Schlüter, M., Tavoni, A., & Levin, S. (2016). Robustness of norm-driven cooperation in the commons. Proceedings of the Royal Society B: Biological Sciences, 283(1822), 20152431] (https://royalsocietypublishing.org/doi/full/10.1098/rspb.2015.2431), an ABM of [Tavoni, A., Schlüter, M., & Levin, S. (2012). The survival of the conformist: social pressure and renewable resource management. Journal of theoretical biology, 299, 152-161.](https://www.sciencedirect.com/science/article/abs/pii/S002251931100347X)

The goal of these simulations is to explore the power and pitfalls of ["link reciprocity"](https://yins.yale.edu/illustrative-projects/social-networks-can-be-used-increase-human-cooperation) in common-pool resource games. I (Elle) am trying to make the argument that though this strategy for network dynamics "works" in games like the Prisoner's Dilemma, it is not particularly helpful in other contexts and may even constitute a social trap. To link this work to the rest of my PhD, we are going to incorperate additional tie rewiring strategies and network mechanisms that are informed by the broader risk hypthesis lit and explore their impact (e.g., in potentially averting the collapse of the system into a defector regime in critical parameter space). 

This builds upon work looking at the impact of static networks in this model, including: 
- Chung, N. N., Chew, L. Y., & Lai, C. H. (2013). Influence of network structure on cooperative dynamics in coupled socio-ecological systems. EPL (Europhysics Letters), 104(2), 28003.
- Sugiarto, H. S., Chung, N. N., Lai, C. H., & Chew, L. Y. (2015). Socioecological regime shifts in the setting of complex social interactions. Physical Review E, 91(6), 062804.
- Sugiarto, H. S., Lansing, J. S., Chung, N. N., Lai, C. H., Cheong, S. A., & Chew, L. Y. (2017). Social cooperation and disharmony in communities mediated through common pool resource exploitation. Physical Review Letters, 118(20), 208301.

The starting point is the lone study (to my knowledge) that incorporated network dynamics: Min, Y., Du, Y., & Jin, C. (2018). The effect of link rewiring on a coevolutionary common pool resource game. Physica A: Statistical Mechanics and its Applications, 512, 935-944. Some claims made in this paper are problematic, for reasons outlined in main_min.R (which may or may not be up to date). For example, and as illustrated below, the finding that "partner switching can stabilize the cooperation when the resource has high inflow or reproductivity" (top; a = <k> = 15 and b = <k> = 90) is driven by the decreasing probability of strategies being updated as the probability of partner switching increases. When the probability of tie rewiring is seperated from the probability of strategy updating (bottom pannel, rewiring probability ~ 1), the effect no longer holds...  

![replicating_min_figure3](https://user-images.githubusercontent.com/48939952/215665702-96e44330-f34d-4dda-99f7-325acc7d92ca.png)

The probability of strategy updating in Min et al. was influenced by this tradeoff, but to a lesser degree, as it is also determined by the normalised payoff difference. This is clearest by comparing the proportion of runs with strategy updating (first series of plots below) and, especially, the proportion of runs with strategy updating *where strategy updating was possible* (second series of plots). 
  
![proportion_runs_strat_updating](https://user-images.githubusercontent.com/48939952/215663568-79150137-4008-4dbe-aea3-8fb25b2f038e.png)

![prob_stratupdating_whenpossible](https://user-images.githubusercontent.com/48939952/215665343-d8117cb0-4312-48f7-b341-5499fc8e297a.png)

[insert summary of dynamics before collapse into defector equilibrium around RI = 50 in <k> = 15 and RI ~53 in <k> = 90 ]   

After this, will (1) introduce different tie updating strategy that aligns with the risk hypothesis and (2) give ties a new purpose (monitoring + sanctioning) 

  
Floyd presented some of this work... 

## SUNBELT 2022
### "Modelling the evolution of cooperation in social-ecological systems: do we need to incorporate network processes or is network structure sufficient?"

Abstract: 
How do cooperative solutions to environmental social dilemmas evolve in populations of rational individuals and what mechanisms maintain them overtime? An important step towards answering this question came with the recognition that ecological dynamics cannot be abstracted away in models of the evolution of cooperation for environmental management, as it is the combination of social and ecological conditions that determine the incentives people have to cooperate, by restraining their natural resource use to sustainable levels, or to defect, by engaging in excessive consumption. A second advancement was made with the move away from assuming that populations are well-mixed when modeling these systems – that is, that they can be represented by complete graphs, where everyone has the opportunity to interact with everyone else and does so with equal probability  – to consider structured populations, including scenarios where agents can rewire their social networks. However, in the pursuit of simple models, some of the work in this vein has neglected the social network and psychological literatures. It remains an open question, then, whether incorporating more realistic network processes will impact the outcomes observed in these models and challenge any of the basic laws they have put forward, some of which are at odds with insights from other disciplines (e.g., that cooperation will evolve if the benefit-to-cost ratio is greater than the average degree and that clustering is a viable solution to many cooperation problems). In this presentation, we begin to address this question with a classic social-ecological systems model, developed by Tavoni, Schlüter, and Levin (2012) and subsequently picked up complex systems scientists (e.g., Sugiarto, Chung, Lai, & Chew, 2015). It incorporates social pressure to comply with cooperative norms regarding water use but, at present, some unrealistic assumptions about the impact of social sanctioning in networks. 


## MARCH 2023 

Sanity checks, requested by Michael, indicate that the parameterisation in Min et al., is problematic for an additional reason. The resource inflow parameter, c, and density dependent depreciation factor, d, need to be the same for the resource dynamics to work as intended. It is unclear whether Min et al. misrepored their method or implemented this incorrectly. Either way, I made a mistake again by following what has been done previously. If this is to be used in my thesis, all of the models need to be re-run (which will take MONTHS) and I'm also going to have to refine the code before doing so. 




