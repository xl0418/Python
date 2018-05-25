import numpy as np
from Trait_sim_in_branches_stat import traitsim
from ABC_MCMC import calibrication,MCMC_ABC


# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
                  theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=1)


# Calibriation step
cal_size = 20000
priorpar = [0.2,0.5,0.1,0.4]
collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, mode='nor')
np.savetxt("/home/p274981/Python_p2/calibration2w_3chains.txt",collection)
#collection = np.loadtxt("/home/p274981/Python_p2/testcal.txt")
