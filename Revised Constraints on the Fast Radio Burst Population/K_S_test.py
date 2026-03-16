from scipy.stats import kstest, ks_2samp
from FRBpopulation.MCMC import Chi_analyze
import FRBpopulation.prepare.CHIME_Lin as CL
from FRBpopulation.MCMC import mcmc
import numpy as np
from scipy.integrate import quad
import FuncZou as fz
import matplotlib.pyplot as plt
from setup import *


def cdff(lgfv):
    return -1.5*np.log10(10**lgfv)

plt.plot(lgfv_space, cdff(lgfv_space))
plt.show()