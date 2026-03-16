from FRBpopulation.setup import *
from FRBpopulation.prepare import dm_model
from FRBpopulation import FuncZou
import matplotlib.pyplot as plt
import pandas as pd
import xlwt

def readsample(i):
    df = pd.read_excel(r'F:\pythonProject1\FRBpopulation\data\selected_samples.xlsx', usecols=[i])
    df_li = df.values.tolist()
    return np.array(df_li)

dm_frb = readsample(0)
z = []
for dm in dm_frb:
    z0 = FuncZou.get_sample(dm_model.p_z_at_dm, 0.001, 3, 1, 100, 1000, dm, sigma_host0, emu0, f_IGM_p, alpha0, F0)
    z.append(z0[0])
print(z)

myworkplace = xlwt.Workbook()
sheet = myworkplace.add_sheet('sheet')
sheet.write(0,0,'z')
for i in range(len(z)):
    sheet.write(i+1,0,z[i])
myworkplace.save(r'F:\pythonProject1\FRBpopulation\data\z_samples.xlsx')
