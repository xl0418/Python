import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.pylab import *


x = np.array(range(2001))
y1 = trait_BH[x,0]
y2 = trait_BH[x,1]

y3 = trait_RI[x,0]
y4 = trait_RI[x,1]

size1 = population_BH[x,0]*10
size2 = population_BH[x,1]*10
size3 = population_RI[x,0]*10
size4 = population_RI[x,1]*10

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
ax01.set_xlim(0, 2001)
ax02.set_xlim(0, 2001)
# ax03.set_xlim(0,5.0)
# ax04.set_xlim(0,5.0)

# Turn on grids
ax01.grid(True)
ax02.grid(True)

ax01.set_xlabel("")
ax01.set_ylabel("Trait Value")
ax02.set_xlabel("Generation")
ax02.set_ylabel("Trait Value")


line1, = ax01.plot([],[], 'b-')
line2, = ax01.plot([],[], 'r-')
line3, = ax02.plot([],[], 'b-')
line4, = ax02.plot([],[], 'r-')

scatter1 = ax01.scatter([],[],  c = 'b', alpha = 0.5)
scatter2 = ax01.scatter([],[], c = 'r', alpha = 0.5)
scatter3 = ax02.scatter([],[],  c = 'b', alpha = 0.5)
scatter4 = ax02.scatter([],[], c = 'r', alpha = 0.5)



def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    scatter1.set_offsets([])
    scatter2.set_offsets([])
    scatter3.set_offsets([])
    scatter4.set_offsets([])
    # scatter1.set_sizes([])
    # scatter2.set_sizes([])
    # scatter3.set_sizes([])
    # scatter4.set_sizes([])
    return line1, line2, line3, line4, scatter1, scatter2, scatter3, scatter4

def animate(i):
    i = (i+1)%(len(x)+1)
    time_text.set_text(time_template % (i))
    data1 = np.hstack((x[i], y1[i]))
    data2 = np.hstack((x[i], y2[i]))
    data3 = np.hstack((x[i], y3[i]))
    data4 = np.hstack((x[i], y4[i]))
    line1.set_data(x[0:i], y1[0:i])
    line2.set_data(x[0:i], y2[0:i])
    line3.set_data(x[0:i], y3[0:i])
    line4.set_data(x[0:i], y4[0:i])
    scatter1.set_offsets(data1)
    scatter2.set_offsets(data2)
    scatter3.set_offsets(data3)
    scatter4.set_offsets(data4)
    scatter1.set_sizes(size1[0:i])
    scatter2.set_sizes(size2[0:i])
    scatter3.set_sizes(size3[0:i])
    scatter4.set_sizes(size4[0:i])
    return line1, line2, line3, line4, scatter1, scatter2, scatter3, scatter4

##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 0.9, '', transform=ax01.transAxes)

ani = animation.FuncAnimation(f0, animate, interval= 1, repeat=False, blit=False, init_func=init)
plt.show()