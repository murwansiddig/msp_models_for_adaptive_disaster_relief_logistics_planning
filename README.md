# Multi-stage Stochastic Programming Methods for Adaptive Disaster Relief Logistics Planning
## authors:
  1. Murwan Siddig: RWTH Aachen University, [siddig@dpo.rwth-aachen.de](siddig@dpo.rwth-aachen.de)
  2. Yongjia Song: Clemson University, [yongjis@clemson.edu](yongjis@clemson.edu)
## abstract:
We consider a logistics planning problem of prepositioning relief items in preparation for an impending hurricane landfall. This problem is modeled as a multiperiod network flow problem where the objective is to minimize the logistics cost of operating the network and the penalty for unsatisfied demand. We assume that the demand for relief items can be derived from the hurricane's predicted intensity and landfall location, which evolves according to a Markov chain. We consider this problem in two different settings, depending on whether the time of landfall is assumed to be deterministic (and known a priori) or random. For the former case, we introduce a fully adaptive MSP model that allows the decision-maker to adjust the prepositioning decisions, sequentially, over multiple stages, as the hurricane's characteristics become clearer. For the latter case, we extend the MSP model with a random number of stages introduced in Guigues (2021), to the case where the underlying stochastic process is assumed to be stage-wise dependent. We benchmark the performance of the MSP models with other approximation policies such as the static and rolling-horizon two-stage stochastic programming approaches. Our numerical results provide key insight into the value of MSP, in disaster relief logistics planning.

## code:
This repository contain two main repositories:
- hurricane_logistics_deterministic: code for situation when the hurricane landfall is *deterministic*  
- hurricane_logistics_random: code for situation when the hurricane landfall is *random*
