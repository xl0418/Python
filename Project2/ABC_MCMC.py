import numpy as np
from Trait_sim_in_branches_stat import traitsim, drawplot, dotplot
import scipy.stats

def single_trait_sim(par):
    sim = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par[0],gamma_K2=par[0],a = par[1],r = 1,theta = 0,K = 5000
                   , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate = 0)
    return sim

def PosNormal(mean, sigma):
    x = np.random.normal(mean,sigma,1)
    return(x if x>=0 else PosNormal(mean,sigma))

def calibrication(samplesize, priorpar, obs, mode = 'uni'):
    collection = np.zeros(shape=(samplesize,4))
    if mode == 'uni':
        uniform_gamma = np.random.uniform(priorpar[0],priorpar[1],samplesize)
        uniform_a = np.random.uniform(priorpar[2],priorpar[3],samplesize)
        do = True
    elif mode == 'nor':
        uniform_gamma = np.zeros(samplesize)
        uniform_a = np.zeros(samplesize)
        for i in range(samplesize):
            uniform_gamma[i] = PosNormal(priorpar[0],priorpar[1])
            uniform_a[i] = PosNormal(priorpar[2],priorpar[3])
        do = True
    elif mode == 'self':
        uniform_gamma = np.repeat(priorpar[0],samplesize)
        uniform_a = np.repeat(priorpar[1],samplesize)
        do = True
    else:
        print('Please indicate one mode!')
        uniform_a = 0
        uniform_gamma = 0
        do = False

    if do == True:
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

def MCMC_ABC(startvalue, iterations,delta,obs,sort,priorpar,mode = 'uni'):
    MCMC = np.zeros(shape=(iterations+1,2))
    MCMC[0,] = startvalue
    par_jump = np.empty(2)
    if mode == 'uni':
        for i in range(iterations):
            par_jump[0] = abs(np.random.normal(loc=MCMC[i,0], scale= 0.01 ))
            par_jump[1] = abs(np.random.normal(loc=MCMC[i,1], scale= 0.01 ))

            if (ABC_acceptance(par_jump,delta = delta, obs = obs,sort = sort)):
                MCMC[i+1,] = par_jump
                print("MCMC : %d Accepted" % (i+1))

            else:
                MCMC[i + 1,] = MCMC[i ,]
                print("MCMC : %d Rejected" % (i+1))
    elif mode == 'nor':
        for i in range(iterations):
            par_jump[0] = abs(np.random.normal(loc=MCMC[i, 0], scale=0.01))
            par_jump[1] = abs(np.random.normal(loc=MCMC[i, 1], scale=0.01))

            pro = np.random.uniform(0,1,1)[0]
            pro_gamma1 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(par_jump[0])
            pro_gamma2 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(MCMC[i ,0])
            pro_a1 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(par_jump[1])
            pro_a2 = scipy.stats.norm(priorpar[0], priorpar[1]).pdf(MCMC[i ,1])

            pro_ratio = (pro_gamma1*pro_a1)/(pro_gamma2*pro_a2)
            accept_creterion = np.min([1,pro_ratio])
            if ABC_acceptance(par = par_jump, delta=delta, obs=obs, sort=sort) and (pro <= accept_creterion):
                MCMC[i + 1,] = par_jump
                print("MCMC : %d Accepted" % (i + 1))
            else:
                MCMC[i + 1,] = MCMC[i,]
                print("MCMC : %d Rejected" % (i + 1))

    return MCMC

