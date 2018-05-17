import matplotlib.animation as animation
from matplotlib.pylab import *
import numpy as np
import matplotlib.pyplot as plt
from JSAnimation import IPython_display



# Function ga: natural selection
def ga(gamma, theta, zi, r):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] =  np.sum(r * np.exp(-gamma * (np.array(theta) - zi[n1]) ** 2)) # np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret
#Derivative of ga
def ga_der(gamma, theta, zi, r):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] =  np.sum(2 * r *gamma * (np.array(theta)-zi[n1]) * np.exp(-gamma * (np.array(theta) - zi[n1]) ** 2)) # np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Function beta: competition part
def beta(a, zi, zj, nj):
    zi_ret = np.ndarray((1,len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0,n1] = np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret
# Function sigma: the derivative of beta with respect to zi
def sigma(a, zi, zj, nj):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = np.sum(2 * a * (zi[n1]-np.array(zj)) * np.exp( -a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret

# Dynamic carrying capacity
def Kd(gamma_K, theta, zi, K):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] = max(np.sum(K * np.exp(-gamma_K * (
                    np.array(theta) - zi[n1]) ** 2)),1)  # np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret
#Derivative of Kd
def Kd_der(gamma_K, theta, zi, K):
    zi_ret = np.ndarray((1, len(zi)))
    for n1 in range(len(zi)):
        zi_ret[0, n1] =  np.sum(2 * K *gamma_K * (np.array(theta)-zi[n1]) * np.exp(-gamma_K * (np.array(theta) - zi[n1]) ** 2)) # np.sum(np.exp(-a * (zi[n1] - np.array(zj)) ** 2) * np.array(nj))
    return zi_ret
  # evolution time: speciation time

speciate_time = 1000
extinction_time = 2000
evo_time = 3000
total_species = 8

trait_BH1 = np.zeros((evo_time + 1, total_species))
population_BH1 = np.zeros((evo_time + 1, total_species))
trait_RI1 = np.zeros((evo_time + 1, total_species))
population_RI1 = np.zeros((evo_time + 1, total_species))

trait_BH2 = np.zeros((evo_time + 1, total_species))
population_BH2 = np.zeros((evo_time + 1, total_species))
trait_RI2 = np.zeros((evo_time + 1, total_species))
population_RI2 = np.zeros((evo_time + 1, total_species))

trait_RI_fin = np.array([trait_RI1,trait_RI2])
trait_BH_fin = np.array([trait_BH1,trait_BH2])
population_BH_fin = np.array([population_BH1,population_BH2])
population_RI_fin = np.array([population_RI1,population_RI2])



randomseed = [12,18]
a_vector = [0.1,0.2]
for loc in np.arange(2):
    np.random.seed(randomseed[loc])
    theta =   np.array([0])   # optimum of natural selection
    gamma = 0.01# intensity of natural selection
    r = 1  # growth rate
    a = a_vector[loc] # intensity of competition
    K = 5000  # carrying capacity
    # parameter controls the type of competition: 1, competition; -1 cooperation.
    m = 1


    # gamma_K = gamma
    # Existing species matrix
    existing_species = np.zeros(shape = (3,total_species))  # np.matrix([[1,1,0], [1,1,1],[1,0,1]])
    existing_species.fill(1)
    # existing_species[0,2]=existing_species[2,1] = 0

    # Variance of random walk of trait evolution
    delta_trait = 0.1
    # build trait evolution and population evolution matrices
    trait_BH = np.zeros((evo_time + 1, total_species))
    population_BH = np.zeros((evo_time + 1, total_species))
    trait_RI = np.zeros((evo_time + 1, total_species))
    population_RI = np.zeros((evo_time + 1, total_species))

    #  initialize condition for species trait and population
    mu_trait, sigma_trait = 0, 10  # mean and standard deviation
    trait_BH[0] = existing_species[0] * np.random.normal(mu_trait, sigma_trait, total_species)
    trait_RI[0] = trait_BH[0]
    print(trait_BH[0])
    print(trait_RI[0])

    mu_pop, sigma_pop = 50, 10  # mean and standard deviation
    population_BH[0] =existing_species[0] * np.random.normal(mu_pop, sigma_pop, total_species)
    population_RI[0] = population_BH[0]
    # print(population_BH[0])


    speciation_event = np.array([1])


    for i in range(evo_time):
        delta_pop = 0.01

        if(i < speciate_time):
            node = 0
        elif(i < extinction_time and i>= speciate_time):
            node = 1
        else:
            node = 2
        # num_species = node + 2
        index_existing_species = np.where(existing_species[node] == 1)[0]
        K_BH = K
        Gamma_BH = ga(gamma=gamma, theta=theta, zi=trait_BH[i,index_existing_species], r=r)
        Gamma_der_BH = ga_der(gamma=gamma, theta=theta, zi=trait_BH[i,index_existing_species], r=r)

        Beta_BH = beta(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
                       nj=population_BH[i,index_existing_species])
        Sigma_BH = sigma(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
                         nj=population_BH[i,index_existing_species])
        trait_BH[i + 1,index_existing_species] = trait_BH[i,index_existing_species] +Gamma_der_BH\
                                                 * ( 1 -  Beta_BH / K_BH) + m * Gamma_BH  * Sigma_BH / K_BH+ \
                                                          np.random.normal(0, delta_trait, len(index_existing_species))
        population_BH[i + 1,index_existing_species] = population_BH[i,index_existing_species] * np.exp(Gamma_BH * (1-Beta_BH/K_BH)
                                          + np.random.normal(0, delta_pop, len(index_existing_species)))
        population_BH[i + 1, np.where(population_BH[i + 1] < 1)] = 0

        Gamma_RI =r
        # Gamma_der_RI = ga_der(gamma=gamma, theta=theta, zi=trait_RI[i,index_existing_species], r=r)
        K_RI = Kd(gamma_K=gamma, theta=theta, zi=trait_RI[i,index_existing_species], K=K)
        Kd_der_RI = Kd_der(gamma_K=gamma, theta=theta, zi=trait_RI[i,index_existing_species], K=K)

        Beta_RI = beta(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
                       nj=population_RI[i,index_existing_species])
        Sigma_RI = sigma(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
                         nj=population_RI[i,index_existing_species])
        trait_RI[i+1,index_existing_species] = trait_RI[i,index_existing_species] +Gamma_RI*Kd_der_RI\
                                               * (Beta_RI/(K_RI**2)) + m * Gamma_RI  * Sigma_RI / K_RI +\
                                               np.random.normal(0, delta_trait, len(index_existing_species))
        population_RI[i+1,index_existing_species] = population_RI[i,index_existing_species] * np.exp(Gamma_RI*(1-Beta_RI/K_RI)
                                                     + np.random.normal(0, delta_pop, len(index_existing_species)))

        population_RI[i+1,np.where(population_RI[i+1]<1)] = 0

        trait_RI_fin[loc] = trait_RI
        trait_BH_fin[loc] = trait_BH
        population_RI_fin[loc]=population_RI
        population_BH_fin[loc]=population_BH

    # trait_BH_fin.append(trait_BH)
    # population_BH_fin.append(population_BH)
    # population_RI_fin.append(population_RI)



trait_RI = np.column_stack((trait_RI_fin[0],trait_RI_fin[1]))
trait_BH = np.column_stack((trait_BH_fin[0],trait_BH_fin[1]))
population_RI = np.column_stack((population_RI_fin[0],population_RI_fin[1]))
population_BH = np.column_stack((population_BH_fin[0],population_BH_fin[1]))

ext_times_RI = []
ext_times_BH = []
ext_spec_index_RI = np.where(population_RI[evo_time,] == 0)[0]
ext_spec_index_BH = np.where(population_BH[evo_time,] == 0)[0]

for j in ext_spec_index_RI:
    ext_time = np.where(population_RI[:, j] == 0)
    ext_times_RI.append(ext_time[0][0])
for j in ext_spec_index_BH:
    ext_time = np.where(population_BH[:, j] == 0)
    ext_times_BH.append(ext_time[0][0])

total_species = len(trait_BH[0,])
x = np.array(range(evo_time))
BH_traits = []
RI_traits = []
BH_sizes = []
RI_sizes = []
for i in np.arange(total_species):
    BH_trait = trait_BH[x,i]
    RI_trait = trait_RI[x,i]
    BH_size = population_BH[x,i]
    RI_size = population_RI[x,i]
    BH_traits.append(BH_trait)
    RI_traits.append(RI_trait)
    BH_sizes.append(BH_size)
    RI_sizes.append(RI_size)



f0 = figure(num = 0, figsize = (12, 8))#, dpi = 100)
f0.suptitle("Trait Evolution", fontsize=12)
ax01 = subplot2grid((2, 1), (0, 0))
ax02 = subplot2grid((2, 1), (1, 0))

for i in theta:
    ax01.axhline(y = i, color='k', linestyle='--',alpha=0.7)
    ax02.axhline(y = i, color='k', linestyle='--',alpha=0.7)



ax01.set_title('RI-dr')
ax02.set_title('RI-dk')

trait_max = max(np.nanmax(trait_RI),np.nanmax(trait_BH))
trait_min = min(np.nanmin(trait_RI),np.nanmin(trait_BH))


# set y-limits
ax01.set_ylim(trait_min-10,trait_max+10)
ax02.set_ylim(trait_min-10,trait_max+10)
# ax03.set_ylim(-0,5)
# ax04.set_ylim(-10,10)

# sex x-limits
ax01.set_xlim(0, evo_time)
ax02.set_xlim(0, evo_time)
# ax03.set_xlim(0,5.0)
# ax04.set_xlim(0,5.0)

# Turn on grids
ax01.grid(True)
ax02.grid(True)

#
ax01.set_xlabel("")
ax01.set_ylabel("Trait Value")
ax02.set_xlabel("Generation")
ax02.set_ylabel("Trait Value")


popu_BH_spec_texts= []
popu_RI_spec_texts= []

text_y1 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)
text_y2 = np.linspace(0.9, 0.9 - (total_species - 1) * 0.05, total_species)

for i in np.arange(total_species):
    popu_BH_spec_text1 = ax01.text(1,text_y1[i], '', transform = ax01.transAxes)
    popu_RI_spec_text1 = ax02.text(1,text_y2[i],'', transform = ax02.transAxes)
    popu_BH_spec_texts.append(popu_BH_spec_text1)
    popu_RI_spec_texts.append(popu_RI_spec_text1)

BH_lines = []
RI_lines = []
BH_scatters = []
RI_scatters = []
for i in  np.arange(total_species):
    BH_line, = ax01.plot([],[], 'b-')
    RI_line, = ax02.plot([],[], 'b-')
    BH_scatter = ax01.scatter([], [], s=0, c='r', alpha=0.3)
    RI_scatter = ax02.scatter([], [], s=0, c='r', alpha=0.3)

    BH_lines.append(BH_line)
    RI_lines.append(RI_line)
    BH_scatters.append(BH_scatter)
    RI_scatters.append(RI_scatter)

def animate(i):
    i = (i+1)%(len(x)+1)
    time_text.set_text(time_template % (i))
    BH_datas = []
    RI_datas = []
    for j in np.arange(total_species):
        BH_data = np.hstack((x[i], BH_traits[j][i]))
        RI_data = np.hstack((x[i], RI_traits[j][i]))
        BH_datas.append(BH_data)
        RI_datas.append(RI_data)

    # for j in np.arange(total_species):
        BH_lines[j].set_data(x[0:i], BH_traits[j][0:i])
        RI_lines[j].set_data(x[0:i], RI_traits[j][0:i])
        BH_scatters[j].set_offsets(BH_datas[j])
        RI_scatters[j].set_offsets(RI_datas[j])
        BH_scatters[j].set_sizes([BH_sizes[j][i]])
        RI_scatters[j].set_sizes([RI_sizes[j][i]])
        # BH_lines[j].set_label("spec %d" % j)
        # Extinct species being labeled by dashed lines
        # if (population_BH[i, j] == 0):
        #     BH_lines[j].set_dashes([2, 2, 2, 2])
        #     BH_lines[j].set_color("red")
        # if (population_RI[i, j] == 0):
        #     RI_lines[j].set_dashes([2, 2, 2, 2])
        #     RI_lines[j].set_color("red")

        # Animating labels
        popu_BH_spec_texts[j].set_text('POS %d = %.1f' % (j+1, population_BH[i,j]))
        popu_RI_spec_texts[j].set_text('POS %d = %.1f' % (j+1, population_RI[i,j]))
        # if (1 not in ext_spec_index_BH):
        # if (i < 1000):
        #     BH_lines[2].set_data([], [])
        #     RI_lines[2].set_data([], [])
        #     BH_lines[1].set_data(x[0:i], BH_traits[1][0:i])
        #     RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
        # elif (i>= 1000 and i< 2000):
        #     BH_lines[2].set_data(x[1000:i], BH_traits[2][1000:i])
        #     RI_lines[2].set_data(x[1000:i], RI_traits[2][1000:i])
        #     BH_lines[1].set_data(x[0:i], BH_traits[1][0:i])
        #     RI_lines[1].set_data(x[0:i], RI_traits[1][0:i])
        # else:
        #     BH_lines[2].set_data(x[1000:i], BH_traits[2][1000:i])
        #     RI_lines[2].set_data(x[1000:i], RI_traits[2][1000:i])
        #     BH_lines[1].set_data(x[0:2000], BH_traits[1][0:2000])
        #     RI_lines[1].set_data(x[0:2000], RI_traits[1][0:2000])


        if (j in ext_spec_index_BH ):
            end_time_BH = ext_times_BH[np.where(j == ext_spec_index_BH)[0][0]]
            if (i>= end_time_BH):
                BH_lines[j].set_data(x[0:end_time_BH], BH_traits[j][0:end_time_BH])

        if (j in ext_spec_index_RI):
            end_time_RI = ext_times_RI[np.where(j == ext_spec_index_RI)[0][0]]
            if (i>= end_time_RI):
                RI_lines[j].set_data(x[0:end_time_RI], RI_traits[j][0:end_time_RI])
    for cy in np.arange(8):
        BH_lines[cy].set_color("green")
        # BH_lines[2].set_color("green")
        RI_lines[cy].set_color("green")

        # RI_lines[2].set_color("green")

    return BH_lines, RI_lines,BH_scatters,RI_scatters

##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 1, '', transform=ax01.transAxes)
#
ani = animation.FuncAnimation(f0, animate, interval= 1, frames= 3000, repeat=False, blit=False) #, init_func=init)
plt.show()
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save('/Users/dudupig/Google 云端硬盘/Python/Project2/S+C.mp4', writer=writer)