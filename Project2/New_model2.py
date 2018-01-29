import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

theta = 0   # optimum of natural selection
gamma = 0.01 # intensity of natural selection
r = 1  # growth rate
a = 0.01  # intensity of competition
K = 3000  # carrying capacity

def ga(gamma, theta, zi, r):
    return r * np.exp(-gamma * (theta - zi) ** 2)


def beta(a, zi, zj, nj):
    zi_ret = np.ndarray((1,len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0,n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret


def sigma(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(2 * a * (zi[n1]-np.array(zj)) * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret


num_time = 4000
num_species = 3
# trait_N1 = np.zeros((num_time+1, num_species))
# population_N1 = np.zeros((num_time+1, num_species))
trait_N2 = np.zeros((num_time+1, num_species))
population_N2 = np.zeros((num_time+1, num_species))

mu_trait, sigma_trait = 0, 1  # mean and standard deviation
# trait_N1[0] = np.random.normal(mu_trait, sigma_trait, num_species)
trait_N2[0] = np.random.normal(mu_trait, sigma_trait, num_species)
# print(trait_N1[0])
mu_pop, sigma_pop = 50, 10  # mean and standard deviation
# population_N1[0] = np.random.normal(mu_pop, sigma_pop, num_species)
population_N2[0] = np.random.normal(mu_pop, sigma_pop, num_species)
ga_vector = np.vectorize(ga)

for i in range(num_time):
    # Gamma_N1 = ga_vector(gamma = gamma, theta = theta, zi = trait_N1[i], r = r)
    # Beta_N1 = beta(a = a, zi = trait_N1[i], zj = trait_N1[i], nj = population_N1[i])
    # Sigma_N1 = sigma(a = a, zi = trait_N1[i], zj = trait_N1[i], nj = population_N1[i])
    # trait_N1[i+1] = trait_N1[i] + 2*gamma * (theta - trait_N1[i]) * (Gamma_N1 -np.exp(Gamma_N1) * Beta_N1/(K+(np.exp(Gamma_N1)-1) * Beta_N1))\
    #                 + (np.exp(Gamma_N1) - 1) * Sigma_N1 / (K+ (np.exp(Gamma_N1) - 1) * Beta_N1)
    # population_N1[i+1] = population_N1[i] * np.exp(Gamma_N1) * K/ (K + (np.exp(Gamma_N1) - 1) * Beta_N1)
    # population_N1[i+1,np.where(population_N1[i+1]<1)] = 0
    # print("Gamma is ", Gamma_N1)
    # print("First part - second part is ", (Gamma_N1 -np.exp(Gamma_N1) * Beta_N1/(K+(np.exp(Gamma_N1)-1) * Beta_N1)))
    # print("Second part is ",(K+ (np.exp(Gamma_N1) - 1) * Beta_N1))

    # print("Sigma is ", Sigma_N1)
    Gamma_N2 = ga_vector(gamma=gamma, theta=theta, zi=trait_N2[i], r=r)
    Beta_N2 = beta(a=a, zi=trait_N2[i], zj=trait_N2[i], nj=population_N2[i])
    Sigma_N2 = sigma(a=a, zi=trait_N2[i], zj=trait_N2[i], nj=population_N2[i])
    trait_N2[i+1] = trait_N2[i] + 2*gamma * (theta - trait_N2[i]) * Gamma_N2 * (1 - Beta_N2/K) + Gamma_N2  * Sigma_N2 / K \
                  - (2*gamma *(theta - trait_N2[i]) * Gamma_N2 * (1 - Beta_N2/K) * np.exp( Gamma_N2 * (1 - Beta_N2/K)) * Beta_N2
                    + (Gamma_N2 * Beta_N2 / K - 1) * Sigma_N2 * np.exp(Gamma_N2 * (1 - Beta_N2/K)) +Sigma_N2)/ \
                    (K+ (np.exp(Gamma_N2 * (1 - Beta_N2/K)) - 1) * Beta_N2)
    population_N2[i+1] = population_N2[i] * np.exp(Gamma_N2*(1-Beta_N2/K)) * K / (K + (np.exp(Gamma_N2 * (1 - Beta_N2/K))
                            -1) * Beta_N2)
    population_N2[i+1,np.where(population_N2[i+1]<1)] = 0
    print(np.exp(Gamma_N2*(1-Beta_N2/K)) * K / (K + (np.exp(Gamma_N2 * (1 - Beta_N2/K))
                            -1) * Beta_N2))
# print(population[:,1])
#
# x = range(0, 1001)
# y = population[:,1]
# p = plt.plot(x, y, "o")

num_plots = num_species

# Have a look at the colormaps here and decide which one you'd like:
# http://matplotlib.org/1.2.1/examples/pylab_examples/show_colormaps.html
from cycler import cycler
colors = [plt.cm.spectral(i) for i in np.linspace(0, 1, num_species)]

colormap = plt.cm.gist_ncar
plt.gca().set_prop_cycle(cycler('color', colors))

# Plot several different functions...
x = np.arange(num_time+1)
labels = []
plt.subplot(2, 2, 1)
for i in range(1, num_plots+1):
    plt.plot(x, trait_N2[:,i-1])

plt.subplot(2, 2, 2)
for i in range(1, num_plots + 1):
    plt.plot(x, population_N2[:,i-1])
plt.subplot(2, 2, 3)
sns.distplot(trait_N2[num_time,:], hist=False, rug=True)
#
# plt.subplot(2, 2, 4)
# sns.distplot(population_N2[num_time,:], hist=False, rug=True);

# I'm basically just demonstrating several different legend options here...
# plt.legend(labels, ncol=4, loc='upper center',
#            bbox_to_anchor=[0.5, 1.1],
#            columnspacing=1.0, labelspacing=0.0,
#            handletextpad=0.0, handlelength=1.5,
#            fancybox=True, shadow=True)
#
plt.show()