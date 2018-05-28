import numpy as np
import matplotlib.pyplot as plt
from Trait_sim_in_branches_stat import traitsim
from ABC_MCMC import calibrication,MCMC_ABC
import pylab as P
import matplotlib.mlab as mlab
from sklearn.neighbors import KernelDensity
import pymc3 as pm
# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=0)


# Calibrication step
cal_size = 10000
# TEST1: Uniform prior distribution example
priorpar = [0.0001,1,0.0001,1]
collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs)
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt")


#TEST2: Normal prior distribution example
priorpar = [0.1,0.2,0.1,0.3]
# collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, mode = 'nor')
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/priorresult/calibration2w.txt")

#TEST3: Normal prior distribution with 3 MCMCs
priorpar = [0.2,0.5,0.1,0.4]
# collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs, mode = 'nor')
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/MCMC3/calibration2w_3chains.txt")



collection = filtered_coll

# distance distribution
P.figure()
dis_data = collection[:,[2,3]]
n, bins, patches = P.hist(dis_data, 15, density=1, histtype='bar',
                            color=['crimson', 'burlywood'],
                            label=['distance', 'sorted distance'])
P.legend()
# plt.show()
# Estimate prior distribution of parameters
# Generate random samples from a mixture of 2 Gaussians
# with modes at 5 and 10
data = np.array(collection[:,1])
data = data.reshape(-1,1)
# Plot the true distribution
x = np.linspace(0, 1, 100)[:, np.newaxis]
norm_vals = mlab.normpdf(x, 0.1, 0.4)
plt.plot(x, norm_vals)
# Plot the data using a normalized histogram
plt.hist(data, 50, density=True)
# Do kernel density estimation
kd = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(data)
# Plot the estimated densty
kd_vals = np.exp(kd.score_samples(x))
plt.plot(x, kd_vals,'--r')
# Show the plots
plt.show()



#
threshold = 0.05
num = threshold*cal_size-1
delta = np.sort(collection[:,3])[int(num)]
mn,idx = min( (collection[i,3],i) for i in range(len(collection[:,3])) )
startvalue_par = collection[idx,:2]

filtered_coll = collection[collection[:,3]<=delta]

# ABC_MCMC step
iterations = 100

posterior = MCMC_ABC(startvalue= startvalue_par, iterations = iterations, delta = delta, obs = obs,sort = 1,
                     priorpar=priorpar, mode = 'nor')
np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/posterior.txt",posterior)
# posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/001result10w/posterior.txt")
pm.autocorrplot(posterior)


# Statistic
posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/Normaldistributionresult/posterior_nor.txt")
# Distribution plots for parameters
gamma_samples = posterior[5000::,0]
a_samples = posterior[5000::, 1]
figdis = plt.figure(figsize=(12, 8))

iterations = 20000
plt.subplot(211)
plt.title(r"""Distribution of $\gamma$ with %d samples""" % iterations)

plt.hist(gamma_samples, histtype='stepfilled',
         color = 'darkred', bins=30, alpha=0.8, density=True)
plt.ylabel('Probability Density')


plt.subplot(212)
plt.title(r"""Distribution of $a$ with %d samples""" % iterations)
plt.hist(a_samples, histtype='stepfilled',
         color = 'darkblue', bins=30, alpha=0.8, density=True)
plt.ylabel('Probability Density')
plt.show()

figdis.savefig('c:/Liang/Googlebox/Research/Project2/python_p2/posterior.png', dpi=figdis.dpi)

# Trace plots
figtra = plt.figure(figsize=(12, 6))


# Plot alpha trace
plt.subplot(211)
plt.title(r'Trace of $\gamma$')
plt.plot(gamma_samples, color = 'darkred')
plt.xlabel('Samples'); plt.ylabel('Parameter')

# Plot beta trace
plt.subplot(212)
plt.title(r'Trace of $a$')
plt.plot(a_samples, color='b')
plt.xlabel('Samples'); plt.ylabel('Parameter')
plt.tight_layout(h_pad=0.8)
plt.show()
figtra.savefig('c:/Liang/Googlebox/Research/Project2/python_p2/posterior_trace.png', dpi=figtra.dpi)
