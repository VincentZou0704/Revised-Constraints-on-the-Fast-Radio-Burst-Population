import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoLocator, MultipleLocator, FixedLocator, AutoMinorLocator

gold = [38.12, 4.99, 0, 3.79, 7.21, 7.37]
full = [109.10, 7.28, 15.10, 3.60, 0.63, 0]

goldx = [2.703, 1, 0, 0.70, 1.44, 1.47]
fullx = [4.1, 1.46, 2.127, 0.65, 0.126,0]

type = ['Gold', 'Full']
color = ['red', 'blue']
yticks = [0,1,2,3,4]
yticks_label = [0,5,10,50,100]

modelx = ['SFH', 'PL', 'CSFH', 'CPL', 'TSE', 'TSRD']
marker_styles = ['o','*','^','D','s','p']
fig, ax = plt.subplots()

plt.axhline(1, color='gray', linestyle='-.', label = r'$\Delta$BIC = 5')
plt.axhline(2, color='gray', linestyle='--', label = r'$\Delta$BIC = 10')

plt.plot(modelx, goldx, marker = 'o', label ='Gold', color = 'red')
plt.plot(modelx, fullx, marker = '^', label ='Full', color = 'blue')


plt.yticks(yticks,yticks_label)
plt.ylabel('$\Delta$BIC')

plt.xlabel('model')
plt.legend()
plt.savefig(r'F:\pythonProject1\figure\M1\BICtest.pdf')
plt.show()