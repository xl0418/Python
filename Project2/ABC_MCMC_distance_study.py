from ABC_MCMC import calibrication
from Trait_sim_in_branches_stat import traitsim
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
# Observation parameters [gamma,a]
par_obs = np.array([0.1,0.1])
# Observation generated
obs = traitsim(h = 1, num_iteration=1,num_species=10,gamma1=par_obs[0],gamma_K2=par_obs[0],a = par_obs[1],r = 1,
               theta = 0,K = 5000 , mean_trait=0,dev_trait=20,mean_pop=50,dev_pop=20, num_time=2000,replicate=1)

# self study
cal_size = 100
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

dis_ori = np.zeros(shape=(6,6))
dis_sort = np.zeros(shape=(6,6))
# Load files and compute the means and variances
for i in range(len(gamma_vec)):
    for j in range(len(a_vec)):
        # collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs,mode='self')
        # str = 'c:/Liang/Googlebox/Research/Project2/python_p2/selfcal%dgamma%da.txt' % (i,j)
        # str = '/home/p274981/Python_p2/selfcal-%dgamma-%da.txt' % (i,j)
        str = 'c:/Liang/Googlebox/Research/Project2/python_p2/priorresult/calibration2w.txt'
        data = np.loadtxt(str)

        dis_ori[i,j] = np.mean(data[:,2])
        dis_sort[i,j] = np.mean(data[:,3])

# collection = data

# Multiple plots for distance distributions
# fig, axes = plt.subplots(nrows=6, ncols=6)
# fig.set_figheight(20)
# fig.set_figwidth(20)
# fig.subplots_adjust(hspace = .001, wspace=.001)

dis_ori = np.array([np.random.normal(i) for i in range(36)])
dis_ori_m = dis_ori.reshape(6,6)
dis_sort = np.array([np.random.normal(i) for i in range(36)])
dis_sort_m = dis_sort.reshape(6,6)

sns.set_style('white')
sns.despine(offset=15, trim=True)
sns.set_context("talk")
fig = plt.figure(figsize=(10, 10))
outer_grid = gridspec.GridSpec(10, 10,wspace=0.0, hspace=0.0, top=0.95, bottom=0.05, left=0.05, right=0.95)
a_axes = plt.subplot(outer_grid[:3,4:])
gamma_axes = plt.subplot(outer_grid[4:,:3])

inner_grid = gridspec.GridSpecFromSubplotSpec(6, 6,
                                              subplot_spec=outer_grid[4:,4:], wspace=0.0, hspace=0.0)

gamma_axes.yaxis.tick_right()
gamma_axes.xaxis.tick_top()


palette1 = sns.color_palette("YlGn_d",n_colors=len(gamma_vec))
palette1.reverse()


for i in range(len(gamma_vec)):
    a_axes.plot(range(len(gamma_vec)),dis_ori_m[i,:],color = palette1[i],lw = 3)
    gamma_axes.plot(dis_ori_m[:,i],range(len(gamma_vec)),color = palette1[i],lw = 3)
    a_axes.set_xticklabels([])
    gamma_axes.set_yticklabels([])



palette2 = sns.color_palette("OrRd_d",n_colors=len(gamma_vec))
palette2.reverse()


for i in range(len(gamma_vec)):
    a_axes.plot(range(len(gamma_vec)),dis_sort_m[i,:],color = palette2[i])
    gamma_axes.plot(dis_sort_m[:,i],range(len(gamma_vec)),color = palette2[i])

gamma_axes.set_ylim(gamma_axes.get_ylim()[::-1])
gamma_axes.set_xlim(gamma_axes.get_xlim()[::-1])

axs = {}

count = 0
# distance distribution
for i in range(len(gamma_vec)):
    for j in range(len(a_vec)):
        axs[count] = plt.subplot(inner_grid[i,j])
        # str = '/home/p274981/Python_p2/selfcal-%dgamma-%da.txt' % (i,j)
        # collection = np.load(str)
        axs[count].set_xlim([0,50])
        axs[count].set_ylim([0,0.45])
        dis_data = collection[:,[2,3]]
        axs[count].hist(dis_data, 15, density=1, histtype='bar',
                                    color=[palette1[1], palette2[1]],
                                    label=['distance', 'sorted distance'])
        # axs[count].axvline(0.1,color = 'k',linestyle = '--')
        if count in [0,6,12,18,24]:
            axs[count].set_xticks([])
            axs[count].set_yticks([])

        elif count in [31,32,33,34]:
            axs[count].set_yticks([])
            axs[count].set_xticks([0,10,20,30,40])
            axs[count].set_xticklabels(['',10,20,30,40])
        elif count in [5,11,17,23,29]:
            axs[count].set_xticks([])
            axs[count].set_yticks([0, 0.2, 0.4])
            axs[count].set_yticklabels(['', 0.2, 0.4])
        elif count == 35:
            axs[count].set_xticks([0, 10, 20, 30, 40])
            axs[count].set_xticklabels(['', 10, 20, 30, 40])
            axs[count].set_yticks([0, 0.2, 0.4])
            axs[count].set_yticklabels(['', 0.2, 0.4])
        elif count == 30:
            axs[count].set_yticks([])
            axs[count].set_xticks([0,10,20,30,40])
            axs[count].set_xticklabels(['',10,20,30,40])

        else:
            axs[count].set_yticks([])
            axs[count].set_xticks([])

        if count in [0,1,2,3,4,5]:
            # subtitle_a = 'a = %f' % a_vec[j]
            axs[count].set_title('a = %s' % a_vec[j])
        if count in [0, 6, 12, 18, 24, 30]:
            axs[count].set_ylabel('$\gamma$ = %s' % gamma_vec[i])

        axs[count].yaxis.tick_right()

        count += 1

# P.legend()
