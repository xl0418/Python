from __future__ import division
import os
import sys
import glob
from JSAnimation import IPython_display
from matplotlib import animation
import pymc3 as pm
import numpy as np
import matplotlib.pyplot as plt
# Generate some data

true_gamma = 0.1
true_a = 0.1

posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/001result/posterior.txt")

# Quickly hacked plotting code
samples = 20000
fig = plt.figure(figsize=(6, 6))
i_width = (true_gamma-.1, true_gamma+.3)
s_width = (true_a-.1, true_a+.3)
samples_width = (0, samples/2)
ax1 = fig.add_subplot(221, xlim=i_width, ylim=samples_width)
ax2 = fig.add_subplot(224, xlim=samples_width, ylim=s_width)
ax3 = fig.add_subplot(223, xlim=i_width, ylim=s_width,
                      xlabel='$gamma$',
                      ylabel='a')
fig.subplots_adjust(wspace=0.0, hspace=0.0)
true_point, =ax3.plot([],[],'r+')
line1, = ax1.plot([], [],'r', lw=1)
line2, = ax2.plot([], [], 'r',lw=1)
line3, = ax3.plot([], [], 'o', lw=2, alpha=.1)
line4, = ax3.plot([], [], 'r',lw=1, alpha=.5)
line5, = ax3.plot([], [], 'k', lw=1)
line6, = ax3.plot([], [], 'k', lw=1)
ax1.set_xticklabels([])
ax2.set_yticklabels([])
#path = plt.scatter([], [])
lines = [line1, line2, line3, line4, line5, line6, true_point]

def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(i):
    true_point.set_data([0.1],[0.1])
    line1.set_data(posterior[:,0][i:0:-1], range(len(posterior[:i,0])))
    line2.set_data(range(len(posterior[:i,1])), posterior[:,1][i:0:-1])
    line3.set_data(posterior[0:i,0], posterior[0:i,1])
    line4.set_data(posterior[0:i,0], posterior[0:i,1])
    intercept = posterior[i,0]
    x = posterior[i,1]
    line5.set_data([intercept, intercept], [x, s_width[1]])
    line6.set_data([intercept, i_width[1]], [x, x])
    return lines

animation.FuncAnimation(fig, animate, init_func=init,
                        frames=samples, interval=1, blit=True)