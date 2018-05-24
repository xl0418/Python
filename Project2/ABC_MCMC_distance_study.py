from ABC_MCMC import calibrication
from Trait_sim_in_branches_stat import traitsim
import numpy as np
import pylab as P

# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=1)

# self study
cal_size = 5000
# Uniform prior distribution example
# priorpar = [0.1,0.1]
gamma_vec = np.array([0,0.001,0.01,0.1,0.5,1])
a_vec = np.array([0,0.001,0.01,0.1,0.5,1])
priorpar = np.zeros(2)
for i in range(len(gamma_vec)):
    for j in range(len(a_vec)):
        priorpar[0]=gamma_vec[i]
        priorpar[1] = a_vec[j]
        collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs,mode='self')
        # str = 'c:/Liang/Googlebox/Research/Project2/python_p2/selfcal%dgamma%da.txt' % (i,j)
        str = '/home/p274981/Python_p2/selfcal-%dgamma-%da.txt' % (i,j)
        np.savetxt(str,collection)



# distance distribution
P.figure()
dis_data = collection[:,[2,3]]
n, bins, patches = P.hist(dis_data, 15, normed=1, histtype='bar',
                            color=['crimson', 'burlywood'],
                            label=['distance', 'sorted distance'])
P.legend()
