import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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
num_vec = np.arange(1,101,1)
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

    # Diversity statistics
    stat_rate_trait_BH[j,:] =trait_BH[num_time,:]
    stat_rate_trait_RI[j, :] = trait_RI[num_time, :]
    stat_rate_popu_BH[j, :] = population_BH[num_time, :]
    stat_rate_popu_RI[j, :] = population_RI[num_time, :]
    j += 1

ext_index_BH = np.where(stat_rate_popu_BH == 0)
ext_index_RI = np.where(stat_rate_popu_RI == 0)
statplot_trait_BH = stat_rate_trait_BH
statplot_trait_BH[ext_index_BH[0],ext_index_BH[1]] = np.nan
statplot_trait_BH_sorted = np.sort(statplot_trait_BH)
statplot_trait_RI = stat_rate_trait_RI
statplot_trait_RI[ext_index_RI[0],ext_index_RI[1]] = np.nan
statplot_trait_RI_sorted = np.sort(statplot_trait_RI)


mask_BH = ~np.isnan(statplot_trait_BH_sorted)
filtered_data_BH = [d[m] for d, m in zip(statplot_trait_BH_sorted.T, mask_BH.T)]
mask_RI = ~np.isnan(statplot_trait_RI_sorted)
filtered_data_RI = [d[m] for d, m in zip(statplot_trait_RI_sorted.T, mask_RI.T)]
merge_BH = np.concatenate( filtered_data_BH, axis=0 )
merge_RI = np.concatenate( filtered_data_RI, axis=0 )



fig = plt.figure(1, figsize=(12, 9))
# plt.xticks(rotation=90)
min_BH = np.amin(merge_BH)-5
max_BH = np.amax(merge_BH)+5
min_RI = np.amin(merge_RI)-5
max_RI = np.amax(merge_RI)+5
global_min = min(min_BH, min_RI)
global_max = max(max_BH, max_RI)

gs = gridspec.GridSpec(1, 4)

# Create an axes instance
ax1 = fig.add_subplot(gs[0,:-2])
bh = ax1.boxplot(filtered_data_BH, 0 , "", patch_artist=True)
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


ri = ax1.boxplot(filtered_data_RI, 0 , "", patch_artist=True)

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

# add legends
ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['BH', 'RI'], loc='upper left')
# add title
ax1.set_title('Trait distribution')
# add x label and y label
ax1.set_xlabel('Species ordered by trait values')
ax1.set_ylabel('Trait values')
ax1.set_ylim(global_min,global_max)

# Major ticks every 20, minor ticks every 5
major_ticks_x = np.arange(0, 101, 20)
minor_ticks_x = np.arange(0, 101, 1)
major_ticks_y = np.arange(-30, 31, 10)
minor_ticks_y = np.arange(-30, 31, 1)

ax1.set_xticks(major_ticks_x)
ax1.set_xticklabels(major_ticks_x)
ax1.set_xticks(minor_ticks_x, minor=True)
ax1.set_yticks(major_ticks_y)
ax1.set_yticks(minor_ticks_y, minor=True)

# And a corresponding grid
# ax.grid(which='both')

# Or if you want different settings for the grids:
ax1.grid(which='minor', alpha=0.2)
ax1.grid(which='major', alpha=0.5)


ax2 = fig.add_subplot(gs[0, 2])
# ax2.spines["top"].set_visible(False)
# ax2.spines["right"].set_visible(False)
# ax2.spines["left"].set_visible(False)
# ax2.get_xaxis().set_ticks([])
# ax2.get_yaxis().set_ticks([])
ax2.set_title('BH model')
# add x label and y label
ax2.set_xlabel('Count')
# ax1.set_ylabel('Trait values')
major_ticks_x2 = np.arange(0, 41, 10)

ax2.set_xticks(major_ticks_x2)
ax2.set_xticklabels(major_ticks_x2)
ax2.set_ylim(global_min,global_max)
# ax2.axes.get_xaxis().set_ticklabels([])
ax2.axes.get_yaxis().set_ticklabels([])
ax2.hist(merge_BH, bins = 100, color="#95d0fc", alpha=0.5, orientation='horizontal')


ax3 = fig.add_subplot(gs[0, 3])
# ax3.spines["top"].set_visible(False)
# ax3.spines["right"].set_visible(False)
# ax3.spines["left"].set_visible(False)
ax3.set_title('RI model')
# add x label and y label
ax3.set_xlabel('Count')
major_ticks_x3 = np.arange(0, 41, 10)

ax3.set_xticks(major_ticks_x3)
ax3.set_xticklabels(major_ticks_x3)
ax3.set_ylim(global_min,global_max)
# ax3.axes.get_xaxis().set_ticklabels([])
ax3.axes.get_yaxis().set_ticklabels([])
ax3.hist(merge_RI, bins = 100, color="#fc5a50",alpha=0.5, orientation='horizontal')

import os
script_dir = os.path.dirname('__file__')
results_dir = os.path.join(script_dir, 'resultes/')
name = "species1btime2qsim1b"
sample_file_name = "%s.pdf" % name

if not os.path.isdir(results_dir):
    os.makedirs(results_dir)

plt.savefig(results_dir + sample_file_name)

plt.close(fig)