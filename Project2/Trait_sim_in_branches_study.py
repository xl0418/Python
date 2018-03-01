import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
from Trait_sim_in_branches_stat import traitsim, drawplot, dotplot

# competition rate vector
# a_vec = np.array([0.01,0.05,0.1,0.5,1])
# natural selection rate vector
# gamma_vec = a_vec

# parameter settings
r = 1
theta = 0
K = 3000
# gamma_K = 0.01
num_time = 2000
num_species = 10
num_time_vec = np.array([2000,10000])
num_species_vec = np.array([1,2,3,5,10])
num_iteration = 100
count1 = 1

gamma1 = 0.01
gamma2 = 0
gamma_K1 = 0
gamma_K2 = 0.01

a = 0.01


# statistics for settings
for num_time in num_time_vec:
    count2 = 1
    for num_species in num_species_vec:
        traitdata = traitsim(num_time = num_time, num_species= num_species, num_iteration= num_iteration,
                             gamma1 = gamma1, gamma2 = gamma2, a = a, r= r,K = K, theta = theta, mean_trait= 0, dev_trait=10, mean_pop= 50,
                             dev_pop= 10, gamma_K1=gamma_K1, gamma_K2 = gamma_K2)
        fig = drawplot(traitdata = traitdata)
        par = (num_species,num_time,num_iteration) #,count1,count2
        # detect the current dir
        script_dir = os.path.dirname('__file__')
        results_dir = os.path.join(script_dir, 'resultes/')
        # file names
        name = "species%d-time%d-sim%d-DRvsDK" % par
        file_name = "%s.pdf" % name
        # if dir doesn't exist, create it
        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)
        # save figs
        plt.savefig(results_dir + file_name)
        # close the windows showing figs
        plt.close(fig)

        name = "species%d-time%d-sim%d-DRvsDK-TV" % par
        file_name = "%s.pdf" % name
        dotfig = dotplot(traitdata=traitdata)
        plt.savefig(results_dir + file_name)
        plt.close(dotfig)
        count2 +=1
    count1 +=1


# dotfig = dotplot(traitdata = traitdata)