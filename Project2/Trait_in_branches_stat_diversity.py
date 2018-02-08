import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os


# Natural selection function
def ga(gamma, theta, zi, r):
    return r * np.exp(-gamma * (theta - zi) ** 2)

# Dynamic carrying capacity
def Kd(gamma_K, theta, zi, K):
    return max(K * np.exp(-gamma_K * (theta - zi) ** 2),1)


# Competition function
def beta(a, zi, zj, nj):
    zi_ret = np.ndarray((1,len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0,n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Derivative of the competition function
def sigma(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(2 * a * (zi[n1]-np.array(zj)) * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Trait simulation function under both Beverton-Holt model and Ricker model
def traitsim(num_time, num_species, num_iteration, gamma, gamma_K, a, r, theta,K , mean_trait, dev_trait, mean_pop, dev_pop):
    j = 0   # initialize the iteration number
    num_vec = np.arange(1,(num_iteration+1),1) # iteration vector
    nrow = len(num_vec)    # Row number of the trait evolution history matrix
    stat_rate_trait_BH = np.empty((nrow,num_species))   # trait evolution matrix under BH
    stat_rate_trait_RI = np.empty((nrow,num_species))  # trait evolution matrix under RI

    stat_rate_popu_BH = np.empty((nrow,num_species))  # population evolution matrix under BH
    stat_rate_popu_RI = np.empty((nrow,num_species)) # population evolution matrix under RI
    # vectorize ga function
    ga_vector = np.vectorize(ga)
    Kd_vector = np.vectorize(Kd)

    # loop for preset iteration
    for loop in num_vec:
        np.random.seed(loop)  # set random seed from 1 to max of iteration
        print(j)
        # Matrices for trait and populations under BH and RI
        trait_BH = np.zeros((num_time+1, num_species))
        population_BH = np.zeros((num_time+1, num_species))
        trait_RI = np.zeros((num_time+1, num_species))
        population_RI = np.zeros((num_time+1, num_species))

        # initialize the input for trait values and populations
        mu_trait, sigma_trait = mean_trait, dev_trait  # mean and standard deviation
        trait_BH[0] = np.random.normal(mu_trait, sigma_trait, num_species)
        trait_RI[0] = trait_BH[0]
        mu_pop, sigma_pop = mean_pop, dev_pop  # mean and standard deviation
        population_BH[0] = np.random.normal(mu_pop, sigma_pop, num_species)
        population_RI[0] = population_BH[0]


        # trait evolution simulation
        for i in range(num_time):
            # BH model
            Gamma_BH = ga_vector(gamma=gamma, theta=theta, zi=trait_BH[i], r=r)
            K_BH = Kd_vector(gamma_K=gamma_K, theta=theta, zi=trait_BH[i], K=K)
            Beta_BH = beta(a=a, zi=trait_BH[i], zj=trait_BH[i], nj=population_BH[i])
            Sigma_BH = sigma(a=a, zi=trait_BH[i], zj=trait_BH[i], nj=population_BH[i])
            trait_BH[i + 1] = trait_BH[i] + 2 * gamma * (theta - trait_BH[i]) * Gamma_BH * (
                    1 - np.exp(Gamma_BH) * Beta_BH / (K_BH + (np.exp(Gamma_BH) - 1) * Beta_BH)) \
                              + (np.exp(Gamma_BH) - 1) * Sigma_BH / (K_BH + (np.exp(Gamma_BH) - 1) * Beta_BH)
            population_BH[i + 1] = population_BH[i] * np.exp(Gamma_BH) * K_BH / (K_BH + (np.exp(Gamma_BH) - 1) * Beta_BH)
            population_BH[i + 1, np.where(population_BH[i + 1] < 1)] = 0

            #RI model
            Gamma_RI = ga_vector(gamma=gamma, theta=theta, zi=trait_RI[i], r=r)
            K_RI = Kd_vector(gamma_K=gamma_K, theta=theta, zi=trait_RI[i], K=K)
            Beta_RI = beta(a=a, zi=trait_RI[i], zj=trait_RI[i], nj=population_RI[i])
            Sigma_RI = sigma(a=a, zi=trait_RI[i], zj=trait_RI[i], nj=population_RI[i])
            trait_RI[i+1] = trait_RI[i] + 2 * (theta - trait_RI[i]) * Gamma_RI * (gamma - r*(gamma - gamma_K)*
                                                                                  Beta_RI/K_RI) + Gamma_RI  * Sigma_RI / K_RI
            population_RI[i+1] = population_RI[i] * np.exp(Gamma_RI*(1-Beta_RI/K_RI))
            population_RI[i+1,np.where(population_RI[i+1]<1)] = 0

        # Diversity statistics
        stat_rate_trait_BH[j,:] =trait_BH[num_time,:]
        stat_rate_trait_RI[j, :] = trait_RI[num_time, :]
        stat_rate_popu_BH[j, :] = population_BH[num_time, :]
        stat_rate_popu_RI[j, :] = population_RI[num_time, :]
        j += 1

    return stat_rate_trait_BH, stat_rate_trait_RI, stat_rate_popu_BH, stat_rate_popu_RI





# boxplots and trait distributions
def drawplot(traitdata):
    # read simulated data
    stat_rate_trait_BH = traitdata[0]
    stat_rate_trait_RI = traitdata[1]
    stat_rate_popu_BH = traitdata[2]
    stat_rate_popu_RI = traitdata[3]

    # find out extinct species and remove the responding trait values
    ext_index_BH = np.where(stat_rate_popu_BH == 0)
    ext_index_RI = np.where(stat_rate_popu_RI == 0)
    statplot_trait_BH = stat_rate_trait_BH
    statplot_trait_BH[ext_index_BH[0],ext_index_BH[1]] = np.nan
    statplot_trait_BH_sorted = np.sort(statplot_trait_BH)
    statplot_trait_RI = stat_rate_trait_RI
    statplot_trait_RI[ext_index_RI[0],ext_index_RI[1]] = np.nan
    statplot_trait_RI_sorted = np.sort(statplot_trait_RI)

    # filter the missing data
    mask_BH = ~np.isnan(statplot_trait_BH_sorted)
    filtered_data_BH = [d[m] for d, m in zip(statplot_trait_BH_sorted.T, mask_BH.T)]
    mask_RI = ~np.isnan(statplot_trait_RI_sorted)
    filtered_data_RI = [d[m] for d, m in zip(statplot_trait_RI_sorted.T, mask_RI.T)]
    # convert the list to an array
    merge_BH = np.concatenate( filtered_data_BH, axis=0 )
    merge_RI = np.concatenate( filtered_data_RI, axis=0 )
    # create a fig
    fig = plt.figure(1, figsize=(12, 9))
    # determine the limits for axises
    min_BH = np.amin(merge_BH)-5
    max_BH = np.amax(merge_BH)+5
    min_RI = np.amin(merge_RI)-5
    max_RI = np.amax(merge_RI)+5
    global_min = min(min_BH, min_RI)
    global_max = max(max_BH, max_RI)
    # multiple plots arrangement
    gs = gridspec.GridSpec(1, 4)

    # Create an axes instance
    ax1 = fig.add_subplot(gs[0,:-2])
    # boxplot BH trait data
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

    # boxplot RI trait data
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

    # Major ticks, minor ticks for axises
    major_ticks_x = np.arange(0, 101, 20)
    minor_ticks_x = np.arange(0, 101, 1)
    major_ticks_y = np.arange(-30, 31, 10)
    minor_ticks_y = np.arange(-30, 31, 1)

    ax1.set_xticks(major_ticks_x)
    ax1.set_xticklabels(major_ticks_x)
    ax1.set_xticks(minor_ticks_x, minor=True)
    ax1.set_yticks(major_ticks_y)
    ax1.set_yticks(minor_ticks_y, minor=True)

    # Or if you want different settings for the grids:
    ax1.grid(which='minor', alpha=0.2)
    ax1.grid(which='major', alpha=0.5)

    # trait distribution plot for BH
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
    # major_ticks_x2 = np.arange(0, 41, 10)
    #
    # ax2.set_xticks(major_ticks_x2)
    # ax2.set_xticklabels(major_ticks_x2)
    ax2.set_ylim(global_min,global_max)
    # ax2.axes.get_xaxis().set_ticklabels([])
    ax2.axes.get_yaxis().set_ticklabels([])
    ax2.hist(merge_BH, bins = 100, color="#95d0fc", alpha=0.5, orientation='horizontal')

    # trait distribution plot for RI
    ax3 = fig.add_subplot(gs[0, 3])
    # ax3.spines["top"].set_visible(False)
    # ax3.spines["right"].set_visible(False)
    # ax3.spines["left"].set_visible(False)
    ax3.set_title('RI model')
    # add x label and y label
    ax3.set_xlabel('Count')
    # major_ticks_x3 = np.arange(0, 41, 10)
    #
    # ax3.set_xticks(major_ticks_x3)
    # ax3.set_xticklabels(major_ticks_x3)
    ax3.set_ylim(global_min,global_max)
    # ax3.axes.get_xaxis().set_ticklabels([])
    ax3.axes.get_yaxis().set_ticklabels([])
    ax3.hist(merge_RI, bins = 100, color="#fc5a50",alpha=0.5, orientation='horizontal')

    plt.show()
    return fig


# trait 1 v.s. trait 2 dotplot
def dotplot(traitdata):
    # read simulated data
    stat_rate_trait_BH = traitdata[0]
    stat_rate_trait_RI = traitdata[1]
    stat_rate_popu_BH = traitdata[2]
    stat_rate_popu_RI = traitdata[3]

    # find out extinct species and remove the responding trait values
    ext_index_BH = np.where(stat_rate_popu_BH == 0)
    ext_index_RI = np.where(stat_rate_popu_RI == 0)
    statplot_trait_BH = stat_rate_trait_BH
    statplot_trait_BH[ext_index_BH[0],ext_index_BH[1]] = np.nan
    # statplot_trait_BH_sorted = np.sort(statplot_trait_BH)
    statplot_trait_RI = stat_rate_trait_RI
    statplot_trait_RI[ext_index_RI[0],ext_index_RI[1]] = np.nan
    # statplot_trait_RI_sorted = np.sort(statplot_trait_RI)
    # filter the missing data
    mask_BH = ~np.isnan(statplot_trait_BH)
    filtered_data_BH = [d[m] for d, m in zip(statplot_trait_BH.T, mask_BH.T)]
    mask_RI = ~np.isnan(statplot_trait_RI)
    filtered_data_RI = [d[m] for d, m in zip(statplot_trait_RI.T, mask_RI.T)]
    # convert the list to an array
    merge_BH = np.concatenate(filtered_data_BH, axis=0)
    merge_RI = np.concatenate(filtered_data_RI, axis=0)

    statplot_popu_BH = stat_rate_popu_BH
    statplot_popu_BH[ext_index_BH[0], ext_index_BH[1]] = np.nan
    # statplot_popu_BH_sorted = np.sort(statplot_popu_BH)
    statplot_popu_RI = stat_rate_popu_RI
    statplot_popu_RI[ext_index_RI[0], ext_index_RI[1]] = np.nan
    # statplot_popu_RI_sorted = np.sort(statplot_popu_RI)

    # filter the missing data
    mask_popu_BH = ~np.isnan(statplot_popu_BH)
    filtered_popu_data_BH = [d[m] for d, m in zip(statplot_popu_BH.T, mask_popu_BH.T)]
    mask_popu_RI = ~np.isnan(statplot_popu_RI)
    filtered_popu_data_RI = [d[m] for d, m in zip(statplot_popu_RI.T, mask_popu_RI.T)]
    # convert the list to an array
    merge_popu_BH = np.concatenate(filtered_popu_data_BH, axis=0)
    merge_popu_RI = np.concatenate(filtered_popu_data_RI, axis=0)


    # create a fig
    fig = plt.figure(1, figsize=(12, 9))
    # determine the limits for axises
    min_BH_t = np.amin(merge_BH)-5
    max_BH_t = np.amax(merge_BH)+5
    min_RI_t = np.amin(merge_RI)-5
    max_RI_t = np.amax(merge_RI)+5
    global_min_t = min(min_BH_t, min_RI_t)
    global_max_t = max(max_BH_t, max_RI_t)

    min_BH_p = np.amin(merge_popu_BH) - 5
    max_BH_p = np.amax(merge_popu_BH) + 5
    min_RI_p = np.amin(merge_popu_RI) - 5
    max_RI_p = np.amax(merge_popu_RI) + 5
    global_min_p = min(min_BH_p, min_RI_p)
    global_max_p = max(max_BH_p, max_RI_p)
    # multiple plots arrangement
    gs = gridspec.GridSpec(1, 2)

    # Create an axes instance
    ax1 = fig.add_subplot(gs[0,0])
    # boxplot BH trait data
    bh = ax1.scatter(merge_popu_BH, merge_BH, s=10, c = "#95d0fc", marker='o', alpha = 0.5)

    # add legends
    # ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['BH', 'RI'], loc='upper left')
    # add title
    ax1.set_title('Trait v.s. Population under BH')
    # add x label and y label
    ax1.set_xlabel('Population size')
    ax1.set_ylabel('Trait values')
    ax1.set_ylim(global_min_t,global_max_t)

    # Create an axes instance
    ax2 = fig.add_subplot(gs[0, 1])
    # boxplot BH trait data
    ri = ax2.scatter(merge_popu_RI, merge_RI, s=10, c="#fc5a50", marker='o', alpha=0.5)

    # add legends
    # ax1.legend([bh["boxes"][0], ri["boxes"][0]], ['BH', 'RI'], loc='upper left')
    # add title
    ax2.set_title('Trait v.s. Population under RI')
    # add x label and y label
    ax2.set_xlabel('Population size')
    ax2.set_ylabel('')
    ax2.set_ylim(global_min_t, global_max_t)

    plt.show()
    return fig






# competition rate vector
a_vec = np.array([0.01,0.05,0.1,0.5,1])
# natural selection rate vector
gamma_vec = a_vec

# parameter settings
r = 1
theta = 0
K = 3000
# gamma_K = 0.01
num_time = 2000
num_species = 100
num_iteration = 100
count1 = 1

import matplotlib.backends.backend_pdf

# statistics for settings
for gamma in gamma_vec:
    gamma_K = gamma
    count2 = 1
    for a in a_vec:
        traitdata = traitsim(num_time = num_time, num_species= num_species, num_iteration= num_iteration,
                             gamma = gamma, a = a, r= r,K = K, theta = theta, mean_trait= 0, dev_trait=10, mean_pop= 50,
                             dev_pop= 10, gamma_K=gamma_K)
        fig = drawplot(traitdata = traitdata)
        par = (num_species,num_time,num_iteration,count1,count2)
        # detect the current dir
        script_dir = os.path.dirname('__file__')
        results_dir = os.path.join(script_dir, 'resultes/')
        # file names
        name = "species%d-time%d-sim%d-com%d-nat%d-DK" % par
        file_name = "%s.pdf" % name
        # if dir doesn't exist, create it
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        # save figs
        plt.savefig(results_dir + file_name)
        # close the windows showing figs
        plt.close(fig)

        name = "species%d-time%d-sim%d-com%d-nat%d-DK-TV" % par
        file_name = "%s.pdf" % name
        dotfig = dotplot(traitdata=traitdata)
        plt.savefig(results_dir + file_name)
        plt.close(dotfig)
        count2 +=1
    count1 +=1


# dotfig = dotplot(traitdata = traitdata)