import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# timelist=pd.read_csv('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timelist.csv', sep='',header=None)

timelist = np.genfromtxt('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timelist.csv', delimiter=',')
timebranch = np.genfromtxt('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timebranch.csv', delimiter=',')
timeend = np.genfromtxt('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\timeend.csv', delimiter=',')
traittable = np.genfromtxt('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\traittable.csv', delimiter=',')
traittable = np.delete(traittable, (0), axis=0)
traittable = np.delete(traittable, (0), axis=1)
ltable = np.genfromtxt('C:\\Liang\\Googlebox\\Python\\Project2\\R-tree_sim\\Ltable.csv', delimiter=',')
ltable = np.delete(ltable, (0), axis=0)
ltable = np.delete(ltable, (0), axis=1)
daughter_index = np.absolute(ltable[:,2])
daughter_index = [int(x) for x in daughter_index]

parent_index = np.absolute(ltable[:,1])
parent_index = [int(x) for x in parent_index]

scalor = 10000
evo_timelist = timelist[:,1]
evo_timelist = np.delete(evo_timelist,0)
evo_timelist = max(evo_timelist)-evo_timelist
evo_timelist = evo_timelist * scalor
evo_timelist = evo_timelist.astype(int)

timebranch = timebranch[:,1]
timebranch = np.delete(timebranch,0)
timebranch = [int(x)-1 for x in timebranch]
timeend = timeend[:,1]
timeend = np.delete(timeend,0)
timeend = [int(x)-1 for x in timeend]

theta = 0   # optimum of natural selection
gamma = 0.01 # intensity of natural selection
gamma1 = gamma #Growth rate's gamma
gamma_K2 = gamma # Carrying capacity's gamma; It can be different from growth rate's gamma.
h=1 # Inheritence
r = 1  # growth rate
a = 0.01  # intensity of competition
K = 5000  # carrying capacity
delta_trait = 0.1 #Variance of random walk of trait evolution
delta_pop = .001 # Variance of random walk of population

# Natural selection function
def ga(gamma, theta, zi, r):
    return r * np.exp(-gamma * (theta - zi) ** 2)

# Dynamic carrying capacity
def Kd(gamma_K, theta, zi, K):
    return max(K * np.exp(-gamma_K * (theta - zi) ** 2),1)


# Competition function
def beta(a, zi, zj, nj):
    zi_ret = np.zeros( len(zi))
    for n1 in range(len(zi)):
        zi_ret[n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Derivative of the competition function
def sigma(a, zi, zj, nj):
    zi_ret = np.zeros( len(zi))
    for n1 in range(len(zi)):
        zi_ret[n1] = np.sum(2 * a * (zi[n1]-np.array(zj)) * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# evolution time: speciation time
evo_time = max(evo_timelist)

speciate_time = evo_timelist[timebranch]
extinct_time = evo_timelist[timeend]
extinct_time[np.where(extinct_time == evo_time)[0]] = -1

extinct_time = np.delete(extinct_time, np.where(extinct_time == evo_time)[0])
total_species = len(speciate_time)



# Initialize trait evolution and population evolution matrices
trait_RI_dr = np.zeros((evo_time + 1, total_species))
population_RI_dr = np.zeros((evo_time + 1, total_species))
trait_RI_dk = np.zeros((evo_time + 1, total_species))
population_RI_dk = np.zeros((evo_time + 1, total_species))

#  initialize condition for species trait and population
mu_trait, sigma_trait = 0, 10  # mean and standard deviation
trait_RI_dr[0, (0,1)] = np.random.normal(mu_trait, sigma_trait, 2)
trait_RI_dk[0] = trait_RI_dr[0]
print(trait_RI_dr[0])
mu_pop, sigma_pop = 500, 10  # mean and standard deviation
population_RI_dr[0, (0,1)] = np.random.normal(mu_pop, sigma_pop, 2)
population_RI_dk[0] = population_RI_dr[0]
print(population_RI_dr[0])

# vectorize function ga
ga_vector = np.vectorize(ga)
Kd_vector = np.vectorize(Kd)
# Existing species matrix
existing_species = traittable



for i in range(evo_time):
    print(i)
    num_event = len(np.where(evo_timelist <= i)[0])
    node = num_event - 2
    index_existing_species = np.where(existing_species[node] == 1)[0]

    # RI dynamic r model
    Gamma_RI_dr = ga_vector(gamma=gamma1, theta=theta, zi=trait_RI_dr[i,index_existing_species], r=r)
    K_RI_dr = K
    Beta_RI_dr = beta(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species],
                      nj=population_RI_dr[i,index_existing_species])
    Sigma_RI_dr = sigma(a=a, zi=trait_RI_dr[i,index_existing_species], zj=trait_RI_dr[i,index_existing_species],
                        nj=population_RI_dr[i,index_existing_species])
    trait_RI_dr[i + 1,index_existing_species] = trait_RI_dr[i,index_existing_species] + h * (2 * gamma1 *
                                        (theta - trait_RI_dr[i,index_existing_species]) * Gamma_RI_dr * (1 -
                                                 Beta_RI_dr / K_RI_dr) + Gamma_RI_dr * Sigma_RI_dr / K_RI_dr
                                               + np.random.normal(0, delta_trait, len(index_existing_species)))
    population_RI_dr[i + 1,index_existing_species] = population_RI_dr[i,index_existing_species] * np.exp(Gamma_RI_dr *
                                          (1 - Beta_RI_dr / K_RI_dr) \
                                                           + np.random.normal(0, delta_pop, len(index_existing_species)))
    # population_RI_dr[i + 1, np.where(population_RI_dr[i + 1] < 1)] = 0
    ext_index_RI_dr = np.where(population_RI_dr[i + 1,index_existing_species] == 0)[0]
    if len(ext_index_RI_dr) > 0:
        trait_RI_dr[i + 1,index_existing_species][ext_index_RI_dr] = 0

    # RI dynamic k model
    Gamma_RI = r
    K_RI = Kd_vector(gamma_K=gamma_K2, theta=theta, zi=trait_RI_dk[i,index_existing_species], K=K)
    Beta_RI = beta(a=a, zi=trait_RI_dk[i,index_existing_species], zj=trait_RI_dk[i,index_existing_species],
                   nj=population_RI_dk[i,index_existing_species])
    Sigma_RI = sigma(a=a, zi=trait_RI_dk[i,index_existing_species], zj=trait_RI_dk[i,index_existing_species],
                     nj=population_RI_dk[i,index_existing_species])
    trait_RI_dk[i + 1,index_existing_species] = trait_RI_dk[i,index_existing_species] + h * (2 * r * gamma_K2 *
                                                    (theta - trait_RI_dk[i,index_existing_species]) / K_RI * Beta_RI
                                               + r * Sigma_RI / K_RI
                                               + np.random.normal(0, delta_trait, len(index_existing_species)))
    population_RI_dk[i + 1,index_existing_species] = population_RI_dk[i,index_existing_species] * np.exp(Gamma_RI * (1 -
                                                                                            Beta_RI / K_RI) \
                                                           + np.random.normal(0, delta_pop, len(index_existing_species)))
    # population_RI_dk[i + 1, np.where(population_RI_dk[i + 1] < 1)] = 0
    ext_index_RI_dk = np.where(population_RI_dk[i + 1,index_existing_species] == 0)[0]
    if len(ext_index_RI_dk) > 0:
        trait_RI_dk[i + 1,index_existing_species][ext_index_RI_dk] = 0
    #
    # Gamma_BH = ga_vector(gamma=gamma, theta=theta, zi=trait_BH[i,index_existing_species], r=r)
    # Beta_BH = beta(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
    #                nj=population_BH[i,index_existing_species])
    # Sigma_BH = sigma(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
    #                  nj=population_BH[i,index_existing_species])
    # trait_BH[i + 1,index_existing_species] = trait_BH[i,index_existing_species] + 2 * gamma * (theta - trait_BH[i,index_existing_species])\
    #                                          * Gamma_BH * ( 1 - np.exp(Gamma_BH) * Beta_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)) +\
    #                   + (np.exp(Gamma_BH) - 1) * Sigma_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH) + \
    #                                          np.random.normal(0, delta_trait, len(index_existing_species))
    # population_BH[i + 1,index_existing_species] = population_BH[i,index_existing_species] * np.exp(Gamma_BH) * K / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)
    # population_BH[i + 1, np.where(population_BH[i + 1] < 1)] = 0
    #
    # Gamma_RI = ga_vector(gamma=gamma, theta=theta, zi=trait_RI[i,index_existing_species], r=r)
    # Beta_RI = beta(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
    #                nj=population_RI[i,index_existing_species])
    # Sigma_RI = sigma(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
    #                  nj=population_RI[i,index_existing_species])
    # trait_RI[i+1,index_existing_species] = trait_RI[i,index_existing_species] + 2*gamma * (theta -
    #                                                                             trait_RI[i,index_existing_species]) \
    #                                        * Gamma_RI * (1 - Beta_RI/K) + Gamma_RI  * Sigma_RI / K + \
    #                                        np.random.normal(0, delta_trait, len(index_existing_species))
    # population_RI[i+1,index_existing_species] = population_RI[i,index_existing_species] * np.exp(Gamma_RI*(1-Beta_RI/K))
    # population_RI[i+1,np.where(population_RI[i+1]<1)] = 0

    if (i+1) in speciate_time:
        spe_event_index = len(np.where(speciate_time <= (i+1))[0])
        trait_RI_dr[i+1,daughter_index[spe_event_index-1]-1] = trait_RI_dr[i+1,parent_index[spe_event_index-1]-1] + \
                                                            np.random.normal(0, 0.01, 1)
        population_RI_dr[i + 1, daughter_index[spe_event_index-1]-1] =1/2 * population_RI_dr[i + 1,
                                                                                       parent_index[spe_event_index-1]-1]
        population_RI_dr[i + 1, parent_index[spe_event_index-1]-1] = 1 / 2 * population_RI_dr[i + 1,
                                                                                        parent_index[spe_event_index-1]-1]
        trait_RI_dk[i + 1, daughter_index[spe_event_index-1]-1] = trait_RI_dk[i + 1,
                                                        parent_index[spe_event_index-1]-1] + np.random.normal(0, 0.01, 1)
        population_RI_dk[i + 1, daughter_index[spe_event_index-1]-1] = 1 / 2 * \
                                                                       population_RI_dk[i + 1, parent_index[spe_event_index-1]-1]
        population_RI_dk[i + 1, parent_index[spe_event_index-1]-1] = 1 / 2 * \
                                                                     population_RI_dk[i + 1, parent_index[spe_event_index-1]-1]

    if (i + 1) in extinct_time:
        extinct_species = int(np.where(extinct_time == (i+1))[0])
        trait_RI_dr[i+1, extinct_species] = None
        population_RI_dr[i+1, extinct_species] = 0
        trait_RI_dk[i + 1, extinct_species] = None
        population_RI_dk[i + 1, extinct_species] = 0


trait_RI_dr[np.where(trait_RI_dr == 0)[0],np.where(trait_RI_dr == 0)[1]] = None
trait_RI_dk[np.where(trait_RI_dk == 0)[0],np.where(trait_RI_dk == 0)[1]] = None




num_plots = total_species

# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
from cycler import cycler
colors = [plt.cm.nipy_spectral(i) for i in np.linspace(0, 1, total_species)]

colormap = plt.cm.gist_ncar
plt.gca().set_prop_cycle(cycler('color', colors))

# Plot several different functions...


trait_dr_tips = trait_RI_dr[evo_time,:][~np.isnan(trait_RI_dr[evo_time,:])]
trait_dk_tips = trait_RI_dk[evo_time,:][~np.isnan(trait_RI_dk[evo_time,:])]
population_RI_dr[population_RI_dr == 0] = None
population_RI_dk[population_RI_dk == 0] = None



x = np.arange(evo_time+1)
labels = []
plt.subplot(2, 4, 1)
for i in range(1, num_plots + 1):
    plt.plot(x, trait_RI_dr[:,i-1])

plt.subplot(2, 4, 2)
for i in range(1, num_plots + 1):
    plt.plot(x, population_RI_dr[:,i-1])
plt.subplot(2, 4, 3)
sns.distplot(trait_dr_tips, hist=False, rug=True)

plt.subplot(2, 4, 4)
sns.distplot(population_RI_dr[evo_time,:], hist=False, rug=True)

plt.subplot(2, 4, 5)
for i in range(1, num_plots + 1):
    plt.plot(x, trait_RI_dk[:,i-1],linewidth=0.1)

plt.subplot(2, 4, 6)
for i in range(1, num_plots + 1):
    plt.plot(x, population_RI_dk[:,i-1])
plt.subplot(2, 4, 7)
sns.distplot(trait_dk_tips, hist=False, rug=True)

plt.subplot(2, 4, 8)
sns.distplot(population_RI_dk[evo_time,:], hist=False, rug=True)

# I'm basically just demonstrating several different legend options here...
# plt.legend(labels, ncol=4, loc='upper center',
#            bbox_to_anchor=[0.5, 1.1],
#            columnspacing=1.0, labelspacing=0.0,
#            handletextpad=0.0, handlelength=1.5,
#            fancybox=True, shadow=True)
#
plt.show()
print(trait_RI_dr[evo_time,:])
#
# histogram=plt.figure()
#
# x =trait_BH[num_time,:]
# y = trait_RI[num_time,:]
# bins = np.linspace(-30, 30, 10)
#
# plt.hist(x, bins, alpha=0.5, label= 'BH')
# plt.hist(y, bins, alpha=0.5, label= 'RI')
# plt.show()
