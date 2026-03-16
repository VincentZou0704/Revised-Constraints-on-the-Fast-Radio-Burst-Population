import matplotlib.pyplot as plt
import numpy as np
from FRBpopulation.setup import *
import pandas as pd


def readCHIME(i):
    df = pd.read_excel(r'F:\pythonProject1\FRBpopulation\data\CHIME.xlsx', usecols=[i])
    df_li = df.values.tolist()
    return df_li

def data_select(gold = True):
    fv_ = np.array(readCHIME('Fluence [Jy ms]')).flatten()
    z_ = np.array(readCHIME('redshift')).flatten()
    lgE_ = np.array(readCHIME('log(E/erg]')).flatten()
    flag_ = np.array(readCHIME('flag')).flatten()
    tag = []

    if gold:
        for i in range(len(fv_)):
            if fv_[i] >= 10**lgfv_min and not np.isnan(z_[i]) and flag_[i] == 1:
               tag.append(i)
    else:
        for i in range(len(fv_)):
            # delete 3 z greater than 3
            if fv_[i] >= 10**lgfv_min and not np.isnan(z_[i]) and z_[i] < 3:
               tag.append(i)

    fv, z, lgE = [], [], []
    for i in range(len(tag)):
        fv.append(fv_[tag[i]])
        z.append(z_[tag[i]])
        lgE.append(lgE_[tag[i]])
    return fv, z, lgE


lgfv0, z0, lgE0 = np.log10(np.array(data_select()[0])),\
                              np.array(data_select()[1]), \
                              np.array(data_select()[2])

lgfv_all, z_all, lgE_all = np.log10(np.array(data_select(gold=False)[0])),\
                              np.array(data_select(gold=False)[1]), \
                              np.array(data_select(gold=False)[2])

if __name__ == '__main__':

    bin = 30
    print('number of FRB:', len(lgfv0), len(lgfv_all))
    print('fv_min: ', np.min(lgfv0), '\t', 'fv_max: ', np.max(lgfv0), '\n',
          'z_min: ', np.min(z0), '\t', 'z_max: ', np.max(z0), '\n',
          'lgE_min: ', np.min(lgE0), '\t', 'lgE_max: ', np.max(lgE0))


    plt.figure()
    plt.hist(lgfv_all, bins = 25, color = '#99CCFF', alpha = 0.55,
             edgecolor = 'white', label='Full sample', range=(lgfv_min, lgfv_max))
    plt.hist((lgfv0), bins = 25, color = '#E11F31', alpha = 0.55,
             edgecolor = 'white', label='Gold sample', range=(lgfv_min, lgfv_max))
    plt.xlabel('log$F_v$')
    plt.ylabel('N(log$F_v$)')
    plt.legend()
    plt.savefig(r'F:\pythonProject1\figure\M1\hist_Fv.pdf')

    plt.figure()
    plt.hist(z_all, bins = bin, color = '#99CCFF', alpha = 0.55,
             edgecolor = 'white', label='Full sample', range=(0, 3))
    plt.hist(z0, bins = bin, color = '#E11F31', alpha = 0.55,
             edgecolor = 'white', label='Gold sample', range=(0, 3))
    plt.xlabel('$z$')
    plt.ylabel('N($z$)')
    plt.legend()
    plt.savefig(r'F:\pythonProject1\figure\M1\hist_z.pdf')

    plt.figure()
    plt.hist(lgE_all, bins=bin, color='#99CCFF', alpha=0.55,
             edgecolor='white', label='Full sample', range=(lgE_min, lgE_max))
    plt.hist(lgE0, bins=bin, color='#E11F31', alpha=0.55,
             edgecolor='white', label='Gold sample', range=(lgE_min, lgE_max))
    plt.xlabel('log$E$')
    plt.ylabel('N(log$E$)')
    plt.legend()
    plt.savefig(r'F:\pythonProject1\figure\M1\hist_E.pdf')

    plt.show()