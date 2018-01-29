import matplotlib.animation as animation
from matplotlib.pylab import *
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(29)

theta = 0   # optimum of natural selection
gamma = 0.01 # intensity of natural selection
r = 1  # growth rate
a = 0.01 # intensity of competition
K = 3000  # carrying capacity

# Function ga: natural selection
def ga(gamma, theta, zi, r):
    return r * np.exp(-gamma * (theta - zi) ** 2)

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
# evolution time: speciation time
speciate_time = 2000
evo_time = 4000
total_species = 3

# build trait evolution and population evolution matrices
trait_BH = np.zeros((evo_time + 1, total_species))
population_BH = np.zeros((evo_time + 1, total_species))
trait_RI = np.zeros((evo_time + 1, total_species))
population_RI = np.zeros((evo_time + 1, total_species))

#  initialize condition for species trait and population
mu_trait, sigma_trait = 0, 5  # mean and standard deviation
trait_BH[0, (0,1)] = np.random.normal(mu_trait, sigma_trait, 2)
trait_RI[0] = trait_BH[0]
print(trait_BH[0])
mu_pop, sigma_pop = 50, 10  # mean and standard deviation
population_BH[0, (0,1)] = np.random.normal(mu_pop, sigma_pop, 2)
population_RI[0] = population_BH[0]
print(population_BH[0])

# vectorize function ga
ga_vector = np.vectorize(ga)

# Existing species matrix
existing_species = np.matrix([[1,1,0], [1,1,1]])
speciation_event = np.array([1])


for i in range(evo_time):
    node = 0 if i < speciate_time else 1
    num_species = node + 2
    index_existing_species = np.where(existing_species[node] == 1)[1]

    Gamma_BH = ga_vector(gamma=gamma, theta=theta, zi=trait_BH[i,index_existing_species], r=r)
    Beta_BH = beta(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
                   nj=population_BH[i,index_existing_species])
    Sigma_BH = sigma(a=a, zi=trait_BH[i,index_existing_species], zj=trait_BH[i,index_existing_species],
                     nj=population_BH[i,index_existing_species])
    trait_BH[i + 1,index_existing_species] = trait_BH[i,index_existing_species] + 2 * gamma * (theta - trait_BH[i,index_existing_species])\
                                             * Gamma_BH * ( 1 - np.exp(Gamma_BH) * Beta_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)) \
                      + (np.exp(Gamma_BH) - 1) * Sigma_BH / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)
    population_BH[i + 1,index_existing_species] = population_BH[i,index_existing_species] * np.exp(Gamma_BH) * K / (K + (np.exp(Gamma_BH) - 1) * Beta_BH)
    population_BH[i + 1, np.where(population_BH[i + 1] < 1)] = 0

    Gamma_RI = ga_vector(gamma=gamma, theta=theta, zi=trait_RI[i,index_existing_species], r=r)
    Beta_RI = beta(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
                   nj=population_RI[i,index_existing_species])
    Sigma_RI = sigma(a=a, zi=trait_RI[i,index_existing_species], zj=trait_RI[i,index_existing_species],
                     nj=population_RI[i,index_existing_species])
    trait_RI[i+1,index_existing_species] = trait_RI[i,index_existing_species] + 2*gamma * (theta - trait_RI[i,index_existing_species]) \
                                           * Gamma_RI * (1 - Beta_RI/K) + Gamma_RI  * Sigma_RI / K
    population_RI[i+1,index_existing_species] = population_RI[i,index_existing_species] * np.exp(Gamma_RI*(1-Beta_RI/K))
    population_RI[i+1,np.where(population_RI[i+1]<1)] = 0

    if (i+1) == speciate_time:
         trait_BH[i+1,2] = trait_BH[i+1,1] + np.random.normal(0, 0.01, 1)
         population_BH[i + 1, 2] =1/2 * population_BH[i + 1, 1]
         population_BH[i + 1, 1] = 1 / 2 * population_BH[i + 1, 1]
         trait_RI[i + 1, 2] = trait_RI[i + 1, 1] + np.random.normal(0, 0.01, 1)
         population_RI[i + 1, 2] = 1 / 2 * population_RI[i + 1, 1]
         population_RI[i + 1, 1] = 1 / 2 * population_RI[i + 1, 1]



evo_time = evo_time+1
x = np.array(range(evo_time))
BH_trait1 = trait_BH[x,0]
BH_trait2 = trait_BH[x,1]
BH_trait3 = trait_BH[x,2]

RI_trait1 = trait_RI[x,0]
RI_trait2 = trait_RI[x,1]
RI_trait3 = trait_RI[x,2]

BH_size1 = population_BH[x,0]
BH_size2 = population_BH[x,1]
BH_size3 = population_BH[x,2]

RI_size1 = population_RI[x,0]
RI_size2 = population_RI[x,1]
RI_size3 = population_RI[x,2]



f0 = figure(num = 0, figsize = (12, 8))#, dpi = 100)
f0.suptitle("Trait Evolution", fontsize=12)
ax01 = subplot2grid((2, 1), (0, 0))
ax02 = subplot2grid((2, 1), (1, 0))
# ax03 = subplot2grid((2, 2), (1, 0), colspan=2, rowspan=1)
# ax04 = ax03.twinx()

ax01.set_title('Beverton-Holt')
ax02.set_title('Ricker')


# set y-limits
ax01.set_ylim(-20,20)
ax02.set_ylim(-20,20)
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
# ax03.grid(True)
#
ax01.set_xlabel("")
ax01.set_ylabel("Trait Value")
ax02.set_xlabel("Generation")
ax02.set_ylabel("Trait Value")
# ax03.set_xlabel("t")
# ax03.set_ylabel("py")
# ax04.set_ylabel("vy")
# ax = plt.axes(xlim = (0,2001), ylim = (-20,20))
# ax2 = plt.axes(xlim = (0,2001), ylim = (-20,20))

BH_line1, = ax01.plot([],[], 'b-')
BH_line2, = ax01.plot([],[], 'r-')
BH_line3, = ax01.plot([],[], 'g-')

RI_line1, = ax02.plot([],[], 'b-')
RI_line2, = ax02.plot([],[], 'r-')
RI_line3, = ax02.plot([],[], 'g-')

BH_scatter1 = ax01.scatter([],[], s = 0, c = 'b', alpha = 0.5)
BH_scatter2 = ax01.scatter([],[], s = 0, c = 'r', alpha = 0.5)
BH_scatter3 = ax01.scatter([],[], s = 0, c = 'g', alpha = 0.5)


RI_scatter1 = ax02.scatter([],[], s = 0, c = 'b', alpha = 0.5)
RI_scatter2 = ax02.scatter([],[], s = 0, c = 'r', alpha = 0.5)
RI_scatter3 = ax02.scatter([],[], s = 0, c = 'g', alpha = 0.5)
#
# def init():
#     BH_line1.set_data([], [])
#     BH_line2.set_data([], [])
#     BH_line3.set_data([], [])
#     RI_line1.set_data([], [])
#     RI_line2.set_data([], [])
#     RI_line3.set_data([], [])
#     BH_scatter1.set_offsets([])
#     BH_scatter2.set_offsets([])
#     BH_scatter3.set_offsets([])
#     RI_scatter1.set_offsets([])
#     RI_scatter2.set_offsets([])
#     RI_scatter3.set_offsets([])
#     return BH_line1, BH_line2, BH_line3, RI_line1, RI_line2,RI_line3,BH_scatter1,BH_scatter2,BH_scatter3,RI_scatter1 \
# , RI_scatter2, RI_scatter3

def animate(i):
    i = (i+1)%(len(x)+1)
    time_text.set_text(time_template % (i))
    BH_data1 = np.hstack((x[i], BH_trait1[i]))
    BH_data2 = np.hstack((x[i], BH_trait2[i]))
    BH_data3 = np.hstack((x[i], BH_trait3[i]))
    RI_data1 = np.hstack((x[i], RI_trait1[i]))
    RI_data2 = np.hstack((x[i], RI_trait2[i]))
    RI_data3 = np.hstack((x[i], RI_trait3[i]))
    BH_line1.set_data(x[0:i], BH_trait1[0:i])
    BH_line2.set_data(x[0:i], BH_trait2[0:i])
    RI_line1.set_data(x[0:i], RI_trait1[0:i])
    RI_line2.set_data(x[0:i], RI_trait2[0:i])
    BH_scatter1.set_offsets(BH_data1)
    BH_scatter2.set_offsets(BH_data2)
    RI_scatter1.set_offsets(RI_data1)
    RI_scatter2.set_offsets(RI_data2)
    BH_scatter3.set_offsets(BH_data3)
    RI_scatter3.set_offsets(RI_data3)
    BH_scatter1.set_sizes([BH_size1[i]])
    BH_scatter2.set_sizes([BH_size2[i]])
    RI_scatter1.set_sizes([RI_size1[i]])
    RI_scatter2.set_sizes([RI_size2[i]])
    BH_scatter3.set_sizes([BH_size3[i]])
    RI_scatter3.set_sizes([RI_size3[i]])
    if i < 2000:
        BH_line3.set_data([], [])
        RI_line3.set_data([], [])

    else:
        BH_line3.set_data(x[2000:i], BH_trait3[2000:i])
        RI_line3.set_data(x[2000:i], RI_trait3[2000:i])

    return BH_line1, BH_line2, BH_line3, RI_line1, RI_line2,RI_line3,BH_scatter1,BH_scatter2,BH_scatter3,RI_scatter1 \
, RI_scatter2, RI_scatter3

##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 0.9, '', transform=ax01.transAxes)

ani = animation.FuncAnimation(f0, animate, interval= 1, frames= 4000, repeat=False, blit=False) #, init_func=init)
plt.show()
#
# Writer = animation.writers['ffmpeg']
# writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save('C:\\Liang\\Googlebox\\Python\\Project2\\S+C.mp4', writer=writer)