import numpy as np
from Trait_sim_in_branches_stat import traitsim, drawplot, dotplot


def single_trait_sim(par):
    sim = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par[0],gamma_K2=par[0],a = par[1],r = 1,theta = 0,K = 5000
                   , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate = 0)
    return sim

def calibrication(samplesize, priorpar, obs):
    collection = np.zeros(shape=(samplesize,4))
    uniform_gamma = np.random.uniform(priorpar[0],priorpar[1],samplesize)
    uniform_a = np.random.uniform(priorpar[2],priorpar[3],samplesize)

    for i in range(samplesize):
        print(i)
        par_cal = np.zeros(2)
        par_cal[0] = uniform_gamma[i]
        par_cal[1] = uniform_a[i]
        sample_cal =  single_trait_sim(par_cal)
        diff =  np.linalg.norm(sample_cal[0]-obs[0])
        diff_sort = np.linalg.norm(np.sort(sample_cal[0])-np.sort(obs[0]))
        collection[i] = np.concatenate((par_cal,[diff],[diff_sort]))
    return collection




def ABC_acceptance(par,delta,obs,sort):
    sample = single_trait_sim(par)
    if sort == 0:
        diff = np.linalg.norm(sample[0] - obs[0])
        if diff<delta:
            return True
        else:
            return False
    else:
        diff_sort = np.linalg.norm(np.sort(sample[0]) - np.sort(obs[0]))
        if diff_sort < delta:
            return True
        else:
            return False

def MCMC_ABC(startvalue, iterations,delta,obs,sort):
    MCMCchain = np.zeros(shape=(iterations+1,2))
    MCMCchain[0,] = startvalue
    par_jump = np.empty(2)
    for i in range(iterations):
        par_jump[0] = abs(np.random.normal(loc=MCMCchain[i,0], scale= 0.01 ))
        par_jump[1] = abs(np.random.normal(loc=MCMCchain[i,1], scale= 0.01 ))
        if (ABC_acceptance(par_jump,delta = delta, obs = obs,sort = sort)):
            MCMCchain[i+1,] = par_jump
            print("MCMC chain: %d Accepted" % (i+1))
        else:
            MCMCchain[i + 1,] = MCMCchain[i ,]
            print("MCMC chain: %d Rejected" % (i+1))

    return MCMCchain

