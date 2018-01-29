import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation

x = np.array(range(2001))
y1 = trait_BH[x,0]
y2 = trait_BH[x,1]

fig = plt.figure()
ax = plt.axes(xlim = (0,2001), ylim = (-30,30))
line1, = ax.plot([],[], 'bo-')
line2, = ax.plot([],[], 'ro-')


def init():
    line1.set_data([], [])
    line2.set_data([], [])

    return line1, line2

def animate(i):
    i = (i+1)%(len(x)+1)
    line1.set_data(x[0:i], y1[0:i])
    line2.set_data(x[0:i], y2[0:i])
    return line1, line2,

ani = animation.FuncAnimation(fig, animate, interval= 1, repeat=True, blit=True, init_func=init)
plt.show()