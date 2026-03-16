import corner
import numpy as np
import pandas as pd
from FRBpopulation.setup import *
import matplotlib.pyplot as plt
from pylab import mpl
import FRBpopulation.FuncZou as fz
import xlwt


Select = False

def getsamples(i):

    ######
    df = pd.read_excel(r'F:\pythonProject1\process\mcmc\CSFH_full.xlsx', usecols=[i], header=None)

    df_li = np.array(df.values.tolist()).flatten()
    return df_li


# labels = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^max$', r'n','$\gamma$', '$z_c$','$\gamma_1$', '$\gamma_2$', '$s$']
labels = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$','n', '$z_c$']



# select
def select(i, rangei):
    index = []
    flat_samples = []
    sample = getsamples(i)
    for j in range(len(sample)):
        if sample[j] >= rangei[0] and sample[j] <= rangei[1]:
            index.append(j)
    for j in range(len(labels)):
        samples = fz.save_by_index(getsamples(j),index)
        flat_samples.append(samples)
    return flat_samples

if Select:
    flat_samples = select(1, (0.,43.2))
else:
    flat_samples = []
    for i in range(len(labels)):
        flat_samples.append(getsamples(i))


flat_samples = np.array(flat_samples).transpose()
print(flat_samples.shape)

def save_mcmc_result(flat_samples, path):
    my_workbook = xlwt.Workbook()
    sheet = my_workbook.add_sheet('mcmc_result')
    for i in range(flat_samples.shape[0]):
        for j in range(flat_samples.shape[1]):
            sheet.write(i, j, flat_samples[i][j])
    my_workbook.save(path)

# save_mcmc_result(flat_samples, r'F:\pythonProject1\process\mcmc\TSRD_full_revised.xlsx')


for i in range(len(labels)):
    perc = np.round(np.percentile(flat_samples[:,i], [15.87, 50, 84.13]),2)
    q = np.round(np.diff(perc), 2)
    print(labels[i], '\t' ,perc[1], '+', q[1], '-', q[0])


# range_x = [1,1,(0.4,0.48),1,1,1,1]


plt.figure()
# mpl.rcParams['font.size'] = 17  # x_ticks
fig = corner.corner(
    flat_samples,
    quantiles=[0.1587, 0.5, 0.8413],
    # range=range_x,
    # label_kwargs={'fontsize': 20},  # axis label
    labels=labels,
    smooth=1, smooth1d=1,
    show_titles = True
)

# font = {'size': 64, 'family' : 'Times New Roman'} # title
# plt.title(r"Gold sample" + '\n' + r'TSRD model', font, x = 0, y = 6)

# plt.savefig(r"F:\pythonProject1\figure\M1\TSRD_gold.pdf")

plt.show()
