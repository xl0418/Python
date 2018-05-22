import elfi
import time
import numpy as np
import scipy.stats
import graphviz
import matplotlib.pyplot as plt
import logging
from Trait_sim_in_branches_stat import traitsim, drawplot, dotplot


par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=0)

# a node is defined by giving a distribution from scipy.stats together with any arguments (here 0 and 2)
gamma = elfi.Prior(scipy.stats.uniform, 0, 1)

# ELFI also supports giving the scipy.stats distributions as strings
a = elfi.Prior('uniform', 0, 1)

def single_trait_sim(gamma,a):
    sim = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=gamma,gamma_K2=gamma,a = a,r = 1,theta = 0,K = 5000
                   , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate = 0)
    return sim

Y = elfi.Simulator(single_trait_sim,gamma,a,observed=obs)

def summary1(x):
    trait = x[0]
    return trait

simdata = elfi.Summary(summary1, Y)

d =  elfi.Distance('euclidean', simdata, obs[0])

elfi.draw(d)
print(d)

import elfi
from elfi.examples import ma2
model = ma2.get_model()
elfi.draw(model)
