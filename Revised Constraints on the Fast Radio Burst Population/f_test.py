from FRBpopulation import FuncZou as fz
from FRBpopulation.setup import *
from FRBpopulation.MCMC import pzE
from scipy.interpolate import interp2d
from scipy import optimize, integrate
import FRBpopulation.prepare.CHIME_Lin as CL
from scipy.stats import chi2
import FRBpopulation.MCMC.Chi_analyze as ca
import xlwt

import matplotlib.pyplot as plt
import numpy as np

# Sample data for the lines
x = np.linspace(0, 10, 100)
y1 = np.sin(x)
y2 = np.cos(x)

# Create the figure and axis
fig, ax = plt.subplots()

# Plot the first line above
ax.plot(x, y1, label='Line 1 (Above)')

# Create a twin x-axis
ax2 = ax.twiny()

# Plot the second line below
ax2.plot(x, y2, label='Line 2 (Below)', color='orange')

# Set labels for both x-axes
ax.set_xlabel('x-axis')
ax2.set_xlabel('x-axis (Secondary)')

# Set y-axis labels for each line
ax.set_ylabel('Line 1 (Above) - Sin(x)')
ax2.set_ylabel('Line 2 (Below) - Cos(x)')

# Add legends for both lines
ax.legend(loc='upper left')
ax2.legend(loc='lower right')

plt.tight_layout()
plt.show()


