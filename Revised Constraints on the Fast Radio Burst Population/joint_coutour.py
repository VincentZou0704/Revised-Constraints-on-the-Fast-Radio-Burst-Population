import corner
import numpy as np
import pandas as pd
from FRBpopulation.setup import *
import matplotlib.pyplot as plt
import getdist.plots as gdplt
from getdist import MCSamples
from scipy.ndimage import gaussian_filter
from matplotlib.ticker import AutoLocator, MultipleLocator, FixedLocator, AutoMinorLocator

type = ['gold', 'full']

def getsamples(i, path):

    ######
    df = pd.read_excel(path, usecols=[i], header=None)

    df_li = np.array(df.values.tolist()).flatten()
    return df_li


# labels = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^max$', r'n','$\gamma$', '$z_c$','$\gamma_1$', '$\gamma_2$', '$s$']
labels = [r'$\alpha$', '$\mathrm{log}E_c$', r'log$F_{\nu,th}^{max}$', r'n']

for i in range(len(type)):
    path = r'F:\pythonProject1\process\mcmc\SFH_' + type[i] + '_revised.xlsx'
    if type[i] == 'gold':
        flat_samples_gold = []
        for j in range(len(labels)):
            flat_samples_gold.append(getsamples(j, path))
    elif type[i] == 'full':
        flat_samples_full = []
        for j in range(len(labels)):
            flat_samples_full.append(getsamples(j, path))


flat_samples_gold = np.array(flat_samples_gold).transpose()
flat_samples_full = np.array(flat_samples_full).transpose()
print(flat_samples_gold.shape, flat_samples_full.shape)


select = True

if not select:
    smooth = np.arange(0, len(labels), 1)
    for i in range(len(smooth)):
        n = int(smooth[i])
        smoothed_samples1 = gaussian_filter(flat_samples_gold[:,n], sigma=0.5)
        smoothed_samples2 = gaussian_filter(flat_samples_full[:,n], sigma=0.5)
        flat_samples_gold[:,n] = smoothed_samples1
        flat_samples_full[:,n] = smoothed_samples2
else:
    smooth = [1,2]
    for i in range(len(smooth)):
        n = int(smooth[i])
        smoothed_samples1 = gaussian_filter(flat_samples_gold[:, n], sigma=2)
        smoothed_samples2 = gaussian_filter(flat_samples_full[:, n], sigma=2)
        flat_samples_gold[:, n] = smoothed_samples1
        flat_samples_full[:, n] = smoothed_samples2


# Create `MCSamples` instances for each chain
samples1 = MCSamples(samples=flat_samples_gold, names= labels)
samples2 = MCSamples(samples=flat_samples_full, names= labels)

# Create a triangle plot with two chains
triangle_plot = gdplt.get_subplot_plotter()
triangle_plot.settings.axes_fontsize = int(2*len(labels)+4) # Set size of ticks
triangle_plot.settings.axes_labelsize = int(2*len(labels)+6)  # Set size of sample labels
triangle_plot.settings.legend_fontsize = int(5*len(labels)+5)  # Set size of legend labels

# must be declared before plot
plt.figure()
triangle_plot.triangle_plot([samples1, samples2], filled=[True, True],
                            # contour_args=[{'smooth': 2.0},{'smooth': 2.0}],
                            contour_colors = ['#FF0000', '#00BFFF'],
                            legend_labels = ['Gold', 'Full'],
                            contour_ls = ['solid', 'dashed'],
                            legend_loc = 'upper right')


# multiplelocator_y = [0.5, 0.05, 2, 5, 2.5, 1]   # TSE
# multiplelocator_y = [0.5, 0.1, 2, 2.5, 2.5, 0.5]   # TSRD
# multiplelocator_y = [0.5, 0.05, 2, 2, 0.25] # CPL
# multiplelocator_y = [0.25, 0.05, 2, 0.5] # PL
# multiplelocator_y = [0.5, 0.025, 2, 0.5]  # CSFH
multiplelocator_y = [0.2, 0.02, 2]  # SFH

for i, row in enumerate(triangle_plot.subplots):

    for j, subplot in enumerate(row):
        if subplot:  # Ensure subplot exists (not None)
            if i == j:
                # major_locator_x = MultipleLocator(multiplelocator[i])
                subplot.xaxis.set_major_locator(AutoLocator())
                subplot.xaxis.set_minor_locator(AutoMinorLocator(4))
            elif i > j:
                major_locator_y = MultipleLocator(multiplelocator_y[i-1])
                subplot.yaxis.set_major_locator(major_locator_y)
                subplot.yaxis.set_minor_locator(AutoMinorLocator(4))
                subplot.xaxis.set_major_locator(AutoLocator())
                subplot.xaxis.set_minor_locator(AutoMinorLocator(4))

# Show the plot
triangle_plot.export(r"F:\pythonProject1\figure\M1\joint_SFH.pdf")  # Export plot to an image file
plt.show()
# plt.savefig(r"F:\pythonProject1\figure\M1\joint_SFH.pdf")
