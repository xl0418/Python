import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

theta = 0   # optimum of natural selection
gamma = 0.01 # intensity of natural selection
r = 1  # growth rate

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


num_time = 2000
j = 0
a_vec = np.arange(0.01, 1.01, 0.01)
num_species = 2
nrow = len(a_vec)
stat_rate = np.empty((nrow,2))
for a in a_vec:
    np.random.seed(29)
    print(j)
    trait_BH = np.zeros((num_time+1, num_species))
    population_BH = np.zeros((num_time+1, num_species))
    trait_RI = np.zeros((num_time+1, num_species))
    population_RI = np.zeros((num_time+1, num_species))

    mu_trait, sigma_trait = 0, 10  # mean and standard deviation
    trait_BH[0] = np.random.normal(mu_trait, sigma_trait, num_species)
    trait_RI[0] = trait_BH[0]
    # print(trait_BH[0])
    mu_pop, sigma_pop = 50, 10  # mean and standard deviation
    population_BH[0] = np.random.normal(mu_pop, sigma_pop, num_species)
    population_RI[0] = population_BH[0]
    ga_vector = np.vectorize(ga)

    for i in range(num_time):
        Gamma_BH = ga_vector(gamma=gamma, theta=theta, zi=trait_BH[i], r=r)
        Beta_BH = beta(a=a, zi=trait_BH[i], zj=trait_BH[i], nj=population_BH[i])
        Sigma_BH = sigma(a=a, zi=trait_BH[i], zj=trait_BH[i], nj=population_BH[i])
        trait_BH[i + 1] = trait_BH[i] + 2 * gamma * (theta - trait_BH[i]) * Gamma_BH * (
                1 - np.exp(Gamma_BH) * Beta_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)) \
                          + (np.exp(Gamma_BH) - 1) * Sigma_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)
        population_BH[i + 1] = population_BH[i] * np.exp(Gamma_BH) * K / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)
        population_BH[i + 1, np.where(population_BH[i + 1] < 1)] = 0
        # print("First part - second part is ", Gamma_BH)
        # print("Second part is ", Gamma_BH * K/ (K + (Gamma_BH - 1) * Beta_BH))
       # print("Sigma is ", Sigma_BH)
        Gamma_RI = ga_vector(gamma=gamma, theta=theta, zi=trait_RI[i], r=r)
        Beta_RI = beta(a=a, zi=trait_RI[i], zj=trait_RI[i], nj=population_RI[i])
        Sigma_RI = sigma(a=a, zi=trait_RI[i], zj=trait_RI[i], nj=population_RI[i])
        trait_RI[i+1] = trait_RI[i] + 2*gamma * (theta - trait_RI[i]) * Gamma_RI * (1 - Beta_RI/K) + Gamma_RI  * Sigma_RI / K
        population_RI[i+1] = population_RI[i] * np.exp(Gamma_RI*(1-Beta_RI/K))
        population_RI[i+1,np.where(population_RI[i+1]<1)] = 0
    # print(population[:,1])
    #
    # x = range(0, 1001)
    # y = population[:,1]
    # p = plt.plot(x, y, "o")

    # Diversity statistics
    stat_rate[j,0] =len(np.where(population_BH[num_time,] == 0)[0])/ num_species
    stat_rate[j, 1] = len(np.where(population_RI[num_time,] == 0)[0]) / num_species
    j += 1



plt.plot(a_vec, stat_rate[:,0], 'r--',a_vec, stat_rate[:,1], 'b--')

