import pandas as pd
import numpy as np
from FRBpopulation.setup import *
import xlwt
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import FRBpopulation.prepare.Fluence as fce


def readCHIME(i):
    df = pd.read_csv(r'F:\pythonProject1\FRBpopulation\data\chimefrbcat1.csv', usecols=[i])
    df_li = df.values.tolist()
    return df_li

# 'bonsai_dm' 'bonsai_snr' 'scat_time' 'dm_exc_ne2001' 'dm_exc_ymw16' 'fluence' 'excluded_flag' 'sub_num=0'

dm_all = np.array(readCHIME('bonsai_dm')) # len = 600
dm_mw_ne2001_all = np.array(readCHIME('dm_exc_ne2001'))
dm_mw_ymw16_all = np.array(readCHIME('dm_exc_ymw16'))
max_mw = np.max([dm_all - dm_mw_ymw16_all, dm_all - dm_mw_ne2001_all], axis=0) + 65

scat_all = np.log10(1000*np.array(readCHIME('scat_time')))
fluence_all = np.log10(np.array(readCHIME('fluence'))+1e-6)     # fluence is scaled by log
excluded_flag = np.array(readCHIME('excluded_flag'))
snr_all = np.array(readCHIME('bonsai_snr'))
repeatornot = readCHIME('repeater_name')
subnum = readCHIME('sub_num')


# tag_test = []
def selection_bias1():
    tag_list1 = []
    for i in range(dm_all.shape[0]):
        # if dm_all[i] > 1.5 * max_mw[i]:
        #     tag_test.append(i)
        if dm_all[i] > 1.5*max_mw[i] and snr_all[i] > 10 and scat_all[i] < 0.8 and \
                fluence_all[i] > 0.5 and excluded_flag[i] == 0 and repeatornot[i][0] == '-9999' and subnum[i][0] == 0:
            tag_list1.append(i)   # to select
    return tag_list1

def selection_bias2():
    tag_list2 = []
    for i in range(dm_all.shape[0]):
        if  fluence_all[i] > -0.5 and repeatornot[i][0] == '-9999' and subnum[i][0] == 0:
            tag_list2.append(i)  # to select
    return tag_list2

tag_list = selection_bias2()
print('number of sample:',len(tag_list))

flence_sample = []
dm_e = []
for i in range(len(tag_list)):
    dm_e.append(dm_mw_ne2001_all[tag_list[i]])
    flence_sample.append(fluence_all[tag_list[i]])

dm_e = np.array(dm_e).transpose()
flence_sample = np.array(flence_sample).transpose()
increase_f = np.sort(flence_sample).reshape(len(tag_list),)
y_f = np.linspace(len(tag_list),1,len(tag_list))

def DM_e_z(z, dm):
    cosmic = coefficient * ob_P * h0_P * f_IGM_p * splz(z) * (1 + alpha0 * z / (1 + z))
    host = 50 / (1 + z)
    return cosmic + host - dm


z_solution = dm_e/855
log_E = fce.iso_E(flence_sample, z_solution)

plt.figure()
plt.plot(increase_f, np.log10(y_f))

plt.figure()
n_count, x_edges, _ = plt.hist(log_E.flatten(), bins= 'auto', alpha = 0)
plt.plot(x_edges[:-1], n_count, linestyle='-', color='red')
plt.show()



def save_DM_z():
    myworkplace = xlwt.Workbook()
    sheet = myworkplace.add_sheet('sheet')
    sheet.write(0,0,'dme')
    sheet.write(0,1,'z')
    sheet.write(0,2,'list_tag')
    for i in range(dm_e.shape[0]):
        sheet.write(i+1,0,dm_e[i])
        sheet.write(i+1,1,z_solution[i])
        sheet.write(i+1,2,tag_list[i])
    myworkplace.save(r'F:\pythonProject1\FRBpopulation\data\selected_samples.xlsx')


# save_DM_z()