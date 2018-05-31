import matplotlib.pyplot as plt
from ABC_MCMC import calibrication,MCMC_ABC
from sklearn.neighbors import KernelDensity
import scipy.stats
import numpy as np

# Calibrication step
cal_size = 20000
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

cal_size = 100000
# TEST4: Uniform prior distribution example
priorpar = [0.0001,1,0.0001,1]
collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs)
# np.savetxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt",collection)
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/testcal.txt")
# collection = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/Calibrtion10w_uni/calibration10w.txt")


# Data filter
threshold = 0.05
num = threshold*cal_size-1
delta = np.sort(collection[:,3])[int(num)]
mn,idx = min( (collection[i,3],i) for i in range(len(collection[:,3])) )
startvalue_par = collection[idx,:2]

filtered_coll = collection[collection[:,3]<=delta]
# Data filter2

delta_unsort = np.sort(collection[:,2])[int(num)]
mn_u,idx_u = min( (collection[i,2],i) for i in range(len(collection[:,2])) )
startvalue_par_u = collection[idx,:2]

filtered2_coll = collection[collection[:,2]<=delta_unsort]



fig, axes = plt.subplots(nrows=3, ncols=2)
fig.set_figheight(10)
fig.set_figwidth(10)
fig.subplots_adjust(hspace = .001, wspace=.001)

axs = axes.ravel()
count = []
data_gamma = np.array(collection[:, 0])
data_a = np.array(collection[:, 1])
data_gamma = data_gamma.reshape(-1, 1)
data_a = data_a.reshape(-1, 1)

dataf_gamma = np.array(filtered_coll[:, 0])
dataf_a = np.array(filtered_coll[:, 1])
dataf_gamma = dataf_gamma.reshape(-1, 1)
dataf_a = dataf_a.reshape(-1, 1)

datafu_gamma = np.array(filtered2_coll[:, 0])
datafu_a = np.array(filtered2_coll[:, 1])
datafu_gamma = datafu_gamma.reshape(-1, 1)
datafu_a = datafu_a.reshape(-1, 1)

x = np.linspace(0, 1, 100)[:, np.newaxis]


# Plot the true distribution
norm_vals_g = scipy.stats.norm.pdf(x, priorpar[0], priorpar[1])
norm_vals_a = scipy.stats.norm.pdf(x, priorpar[2], priorpar[3])
for i in [0,2,4]:
    if i == 0:
        axs[i].plot(x, norm_vals_g)
        axs[i + 1].plot(x, norm_vals_a,label = 'Prior')
    else:
        axs[i].plot(x, norm_vals_g)
        axs[i+1].plot(x, norm_vals_a)

    # Plot the data using a normalized histogram
axs[0].hist(data_gamma, 50, density=True)
axs[2].hist(dataf_gamma, 50, density=True)
axs[4].hist(datafu_gamma, 50, density=True)

axs[1].hist(data_a, 50,label = "Calibration data", density=True)
axs[3].hist(dataf_a, 50, density=True)
axs[5].hist(datafu_a, 50, density=True)


# Do kernel density estimation
kd_g = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(data_gamma)
kdf_g = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(dataf_gamma)
kdfu_g = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(datafu_gamma)

kd_a = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(data_a)
kdf_a = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(dataf_a)
kdfu_a = KernelDensity(kernel='gaussian', bandwidth=0.75).fit(datafu_a)



# Plot the estimated densty
kd_vals_g = np.exp(kd_g.score_samples(x))
kd_vals_gf = np.exp(kdf_g.score_samples(x))
kd_vals_gfu = np.exp(kdfu_g.score_samples(x))

kd_vals_a = np.exp(kd_a.score_samples(x))
kd_vals_af = np.exp(kdf_a.score_samples(x))
kd_vals_afu = np.exp(kdfu_a.score_samples(x))


axs[0].plot(x, kd_vals_g,'--r')
axs[1].plot(x, kd_vals_a,'--r',label = 'Regression')
axs[2].plot(x, kd_vals_gf,'--r')
axs[3].plot(x, kd_vals_af,'--r')
axs[4].plot(x, kd_vals_gfu,'--r')
axs[5].plot(x, kd_vals_afu,'--r')

for i in [1,3,5]:
    axs[i].set_yticks([])
axs[0].set_title('$\gamma$')
axs[1].set_title('a')
axs[0].set_ylabel('Prior-Prior distribution')
axs[2].set_ylabel('Filtered distribution by sorted trait')
axs[4].set_ylabel('Filtered distribution')

axs[1].legend()