from __future__ import division
from matplotlib import animation
import numpy as np
import matplotlib.pyplot as plt

def pretty_array(x):
    return '(%s)' % ', '.join('%.2f' % x_i for x_i in x)

# Load data
true_gamma = 0.1
true_a = 0.1
true_mean = [true_gamma,true_a]
posterior = np.loadtxt("c:/Liang/Googlebox/Research/Project2/python_p2/001result10w/posterior.txt")
N, D = posterior.shape
posterior = posterior[0:N:100,:]

q75_1, q25_1 = np.percentile(posterior[:,0], [75 ,25])
iqr_1 = q75_1 - q25_1

q75_2, q25_2 = np.percentile(posterior[:,1], [75 ,25])
iqr_2 = q75_2 - q25_2

true_iqr = [iqr_1,iqr_2]
# Quickly hacked plotting code
samples = len(posterior[:,0])
fig = plt.figure(figsize=(10, 10))
i_width = (true_gamma-.1, true_gamma+1)
s_width = (true_a-.1, true_a+1)
samples_width = (0, samples/10)
ax1 = fig.add_subplot(221, xlim=i_width, ylim=samples_width)
ax2 = fig.add_subplot(224, xlim=samples_width, ylim=s_width)
ax3 = fig.add_subplot(223, xlim=i_width, ylim=s_width,
                      xlabel='$\gamma$',
                      ylabel='a')
ax1.set_xticks([])
ax2.set_xticks([])
ax1.set_yticks([])
ax2.set_yticks([])

sq_size = .25
j = [0,1]
ax4 = fig.add_axes(
    [1 - .1 - 1.3 * sq_size * (1 - j[0] * D ** -1.) -0.025, 1. - .1 - 1.5 * sq_size * D ** -1 ,
     1.5 * sq_size * D ** -1. - .05, 1.5 * sq_size * D ** -1 - .05])

ax5 = fig.add_axes([0.7625-0.075, 1. - .1 - 1.5 * sq_size * D ** -1 ,
                    1.5 * sq_size * D ** -1. - .05, 1.5 * sq_size * D ** -1. - .05])
text = fig.text(0.55, 0.51,'')

ax4.set_xlabel('$\gamma$')
ax4.set_ylabel('autocorr')
ax4.set_xticks([])
ax4.set_yticks([])
ax4.axis([-15, 15, -.3, 1])

# ax5.acorr(posterior[int(240 / 2.):240:10, 1], detrend=plt.mlab.detrend_mean)
ax5.set_xlabel('a')
ax5.set_xticks([])
ax5.set_yticks([])
ax5.axis([-15, 15, -.3, 1])

fig.subplots_adjust(wspace=0.0, hspace=0.0)

def animate(i):
    ax1.clear()
    ax2.clear()
    ax3.clear()
    ax3.set_xlim(0,1.1)
    ax3.set_ylim(0,1.1)
    ax1.set_xlim(0,1.1)
    ax2.set_ylim(0,1.1)
    ax3.set_xlabel('$\gamma$')
    ax3.set_ylabel(ylabel='a')

    ax3.plot([0.1],[0.1], 'r+')
    ax1.plot(posterior[i::-1,0], range(len(posterior[i::-1,0])), 'r', lw=1)
    ax2.plot(range(len(posterior[i::-1,1])), posterior[i::-1,1], 'r', lw=1)
    ax3.plot(posterior[0:i,0], posterior[0:i,1], 'o', lw=2, alpha=.1)
    ax3.plot(posterior[0:i,0], posterior[0:i,1], 'r', lw=1, alpha=.5)

    intercept = posterior[i, 0]
    x = posterior[i, 1]
    ax3.plot([intercept, intercept], [x, s_width[1]], 'k', lw=1)
    ax3.plot([intercept, i_width[1]],[x, x], 'k', lw=1)
    ax1.set_xticklabels([])
    ax2.set_yticklabels([])

    #autocorr
    if i > 250:
        ax4.clear()
        ax5.clear()

        ax4.set_xlabel('$\gamma$')
        ax4.set_ylabel('autocorr')
        ax4.set_xticks([])
        ax4.set_yticks([])
        ax4.axis([-15, 15, -.3, 1])

        ax5.set_xlabel('a')
        ax5.set_xticks([])
        ax5.set_yticks([])
        ax5.axis([-15, 15, -.3, 1])
        ax4.acorr(posterior[int(i / 2.):i:10, 0], detrend=plt.mlab.detrend_mean)
        ax5.acorr(posterior[int(i / 2.):i:10, 1], detrend=plt.mlab.detrend_mean)

    ## textual information

    # text.set_text('')

    # text.set_text(str)
    str = ''
    if i > 10:

        str += 't = %d\n' % i
        str += 'acceptance rate = %.2f\n\n' % (1. - np.mean(np.diff(posterior[int(i / 2.):i, 0]) == 0.))
        str += 'mean(X) = %s' % pretty_array(posterior[int(i / 2.):i, :].mean(0))
        str += ' / true mean = %s\n' % pretty_array(true_mean)
        # iqr = plt.sort(posterior[int(i / 2.):i, :], axis=0)[int(.25 * (i / 2.), .75 * (i / 2.)), :].T
        q75_1_i, q25_1_i = np.percentile(posterior[int(i / 2.):i, 0], [75, 25])
        # iqr_1_i = q75_1_i - q25_1_i
        q75_2_i, q25_2_i = np.percentile(posterior[int(i / 2.):i, 1], [75, 25])
        # iqr_2_i = q75_2_i - q25_2_i

        # iqr_i = [iqr_1_i, iqr_2_i]

        str += 'IQR($gamma$) = (%.2f, %.2f)' % (q25_1_i, q75_1_i)
        str += ' / true IQR = %.4f\n' % true_iqr[0]
        str += 'IQR(a) = (%.2f, %.2f)' % (q25_2_i, q75_2_i)
        str += ' / true IQR = %.4f\n' % true_iqr[1]

    text.set_text(str)
    plt.draw()


animation.FuncAnimation(fig, animate, frames=(samples-1), interval=1,repeat = False)#, blit=True)