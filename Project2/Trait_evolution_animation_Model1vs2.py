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

line1, = ax01.plot([],[], 'b-')
line2, = ax01.plot([],[], 'r-')
line3, = ax02.plot([],[], 'b-')
line4, = ax02.plot([],[], 'r-')

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    return line1, line2, line3, line4

def animate(i):
    i = (i+1)%(len(x)+1)
    time_text.set_text(time_template % (i))
    line1.set_data(x[0:i], y1[0:i])
    line2.set_data(x[0:i], y2[0:i])
    line3.set_data(x[0:i], y3[0:i])
    line4.set_data(x[0:i], y4[0:i])
    return line1, line2, line3, line4

##
time_template = 'Time = %d G'    # prints running simulation time
time_text = ax01.text(0.05, 0.9, '', transform=ax01.transAxes)

ani = animation.FuncAnimation(f0, animate, interval= 1, repeat=False, blit=False, init_func=init)
plt.show()