import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})
gamma_vec = np.array([0,0.001,0.01,0.1,0.5])
a_vec = np.array([0,0.001,0.01,0.1,0.5,1])
dis_along_gamma_ori = []
dis_along_gamma_sort = []

# Load files and compute the means and variances
# for i in range(len(gamma_vec)):
for i in range(len(gamma_vec)):

    for j in range(len(a_vec)):
        # collection = calibrication(samplesize = cal_size, priorpar = priorpar, obs = obs,mode='self')
        str = 'c:/Liang/Googlebox/Research/Project2/python_p2/Self_dis_study/selfcal-%dgamma-%da.txt' % (i,j)
        # str = '/home/p274981/Python_p2/selfcal-%dgamma-%da.txt' % (i,j)
        # str = 'c:/Liang/Googlebox/Research/Project2/python_p2/priorresult/calibration2w.txt'
        # str = 'd:/Googlebox/Research/Project2/python_p2/priorresult/calibration2w.txt'

        data = np.loadtxt(str)

        dis_along_gamma_ori.append(data[:,2])
        dis_along_gamma_sort.append(data[:,3])
datarow,datacol =data.shape
dis_gamma_ori = np.concatenate( dis_along_gamma_ori, axis=0 )
dis_gamma_sort = np.concatenate( dis_along_gamma_sort, axis=0 )

# collection = data
# gamma_label = np.repeat(['$\gamma$=0','$\gamma$=.001'], 30000)

gamma_label = np.repeat(['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5'], datarow*len(a_vec))
a_label = np.repeat(['a=0','a=.001','a=.01','a=.1','a=.5','a=1'], datarow)
label =np.tile(a_label,len(gamma_vec))
df_ori = pd.DataFrame(dict(dis_gamma_ori=dis_gamma_ori,dis_gamma_sort=dis_gamma_sort,
                           label=label,gamma_label = gamma_label))


#
# # Create the data
# rs = np.random.RandomState(1979)
# x = rs.randn(500)
# g = np.tile(list("ABCDEFGHIJ"), 50)
# df = pd.DataFrame(dict(x=x, g=g))
# m = df.g.map(ord)
# df["x"] += m

# Initialize the FacetGrid object
pal = sns.cubehelix_palette(10, rot=-.25, light=.7)
pal1 = sns.cubehelix_palette(9, rot=-.25, light=.7)
#
# def kde_plot(x, color_pal, **kwargs):
#     cmap = sns.cubehelix_palette(color_pal, rot=-.25, light=.7)
#     sns.kdeplot(x, color_pal=cmap, **kwargs)


g_ori = sns.FacetGrid(df_ori, col = "gamma_label"  ,row="label", hue="label", aspect=15, size=.5,
                      palette=pal,
                      xlim=[0,100])

# Draw the densities in a few steps
g_ori.map(sns.kdeplot, "dis_gamma_ori", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g_ori.map(sns.kdeplot, "dis_gamma_ori", clip_on=False, color="w", lw=2, bw=.2)


g_ori.map(sns.kdeplot, "dis_gamma_sort", clip_on=False, shade=True, alpha=1, lw=1.5, bw=.2)
g_ori.map(sns.kdeplot, "dis_gamma_sort", clip_on=False, color="w", lw=2, bw=.2)
g_ori.map(plt.axhline, y=0, lw=2, clip_on=False)

# # Define and use a simple function to label the plot in axes coordinates
# def label(x, color, label):
#     ax = plt.gca()
#     ax.text(0, .2, label, fontweight="bold", color=color,
#             ha="left", va="center", transform=ax.transAxes)
# #
# # g_ori.axes[0,5].text(0, .2, label, fontweight="bold", color=pal,
# #             ha="left", va="center", transform=g_ori.axes[0,5].transAxes)
# g_ori.map(label, "dis_gamma_sort")

gamma_str = ['$\gamma$=0','$\gamma$=.001','$\gamma$=.01','$\gamma$=.1','$\gamma$=.5']
# gamma_str = ['$\gamma$=0','$\gamma$=.001']

for i in range(len(gamma_str)):
    g_ori.axes[0,i].text(1, 0.8, s = gamma_str[i], fontweight="bold", color = pal[0],
                          ha="right", va="center", transform=g_ori.axes[0,i].transAxes)
a_str = ['a=0','a=.001','a=.01','a=.1','a=.5','a=1']
for i in range(len(a_str)):
    g_ori.axes[i, 4].text(1, .2, s = a_str[i], fontweight="bold", color = pal[i],
                          ha="right", va="center", transform=g_ori.axes[i, 4].transAxes)
# Set the subplots to overlap
g_ori.fig.subplots_adjust(hspace=-.25)

# Remove axes details that don't play will with overlap
g_ori.set_titles("")
g_ori.set(yticks=[])
g_ori.despine(bottom=True, left=True)
# g_ori.set_xlabels("Distance")


