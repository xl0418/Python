import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

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


num_time = 2000
j = 0
num_species = 100
num_vec = np.arange(1,1001,1)
nrow = len(num_vec)
stat_rate_trait_BH = np.empty((nrow,num_species))
stat_rate_trait_RI = np.empty((nrow,num_species))

stat_rate_popu_BH = np.empty((nrow,num_species))
stat_rate_popu_RI = np.empty((nrow,num_species))

for loop in num_vec:
    np.random.seed(loop)
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
    stat_rate_trait_BH[j,:] =trait_BH[num_time,:]
    stat_rate_trait_RI[j, :] = trait_RI[num_time, :]
    stat_rate_popu_BH[j, :] = population_BH[num_time, :]
    stat_rate_popu_RI[j, :] = population_RI[num_time, :]
    j += 1

ext_index_BH = np.where(stat_rate_popu_BH == 0)
ext_index_RI = np.where(stat_rate_popu_RI == 0)
color_BH = dict(boxes='DarkGreen', whiskers='DarkOrange', medians='DarkBlue', caps='Gray')
statplot_trait_BH = stat_rate_trait_BH
statplot_trait_BH[ext_index_BH[0],ext_index_BH[1]] = np.nan
statplot_trait_BH_sorted = np.sort(statplot_trait_BH)
df_BH = pd.DataFrame(statplot_trait_BH_sorted)
# df_BH.plot.box(rot=90, color = color_BH, layout = (2,1), showfliers=False)

color_RI = dict(boxes='DarkRed', whiskers='DarkBlue', medians='DarkBlue', caps='Gray')
statplot_trait_RI = stat_rate_trait_RI
statplot_trait_RI[ext_index_RI[0],ext_index_RI[1]] = np.nan
statplot_trait_RI_sorted = np.sort(statplot_trait_RI)
df_RI = pd.DataFrame(statplot_trait_RI_sorted)
# df_RI.plot.box(rot=90,  color = color_RI, layout = (2,2), showfliers=False)
# Create a figure instance
fig = plt.figure(1, figsize=(9, 6))
# Create an axes instance
ax = fig.add_subplot(111)
bh = ax.boxplot(df_BH.dropna().values, 0 , "", patch_artist=True)
## change outline color, fill color and linewidth of the boxes
for box in bh['boxes']:
    # change outline color
    box.set( color='DarkBlue', linewidth=0.5)
    # change fill color
    box.set( facecolor = '#95d0fc' ,  alpha=0.5)

## change color and linewidth of the whiskers
for whisker in bh['whiskers']:
    whisker.set(color='#95d0fc', linewidth=0.5)

## change color and linewidth of the caps
for cap in bh['caps']:
    cap.set(color='#95d0fc', linewidth=0.5)

## change color and linewidth of the medians
for median in bh['medians']:
    median.set(color='#95d0fc', linewidth=0.5)

## change the style of fliers and their fill
for flier in bh['fliers']:
    flier.set(marker='o', color='#95d0fc', alpha=0.5)


ri = ax.boxplot(df_RI.dropna().values, 0 , "", patch_artist=True)
## change outline color, fill color and linewidth of the boxes
for box in ri['boxes']:
    # change outline color
    box.set( color='#ff000d', linewidth=0.5)
    # change fill color
    box.set( facecolor = '#fc5a50' ,  alpha=0.5)

## change color and linewidth of the whiskers
for whisker in ri['whiskers']:
    whisker.set(color='#fc5a50', linewidth=0.5)

## change color and linewidth of the caps
for cap in ri['caps']:
    cap.set(color='#fc5a50', linewidth=0.5)

## change color and linewidth of the medians
for median in ri['medians']:
    median.set(color='#fc5a50', linewidth=0.5)

## change the style of fliers and their fill
for flier in ri['fliers']:
    flier.set(marker='o', color='#fc5a50', alpha=0.5)
