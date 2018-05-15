import numpy as np
from Trait_sim_in_branches_stat import traitsim
from ABC_MCMC import single_trait_sim,calibrication,ABC_acceptance,MCMC_ABC


# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
                  theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=0)


# Calibrication step
cal_size = 10000
# priorpar = [0.0001,1,0.0001,1]
# collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs)
# np.savetxt("/home/p274981/Python_p2/testcal.txt",collection)
collection = np.loadtxt("/home/p274981/Python_p2/testcal.txt")

threshold = 0.05
num = threshold*cal_size-1
delta = np.sort(collection[:,3])[int(num)]
mn,idx = min( (collection[i,3],i) for i in range(len(collection[:,3])) )
startvalue_par = collection[idx,:2]

# ABC_MCMC step
iterations = 20000

posterior = MCMC_ABC(startvalue= startvalue_par, iterations = iterations, delta = delta, obs = obs,sort = 1)
np.savetxt("/home/p274981/Python_p2/posterior.txt",posterior)
# posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/posterior.txt")

