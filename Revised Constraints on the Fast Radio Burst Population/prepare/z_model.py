import numpy as np
from FRBpopulation import FuncZou as fz
import matplotlib.pyplot as plt
from FRBpopulation.setup import *


def pz(z, model = None, args = None):
    # distribution function pz --- normalisation is not given

    SFH = (1 + z) ** 2.6 / (1 + ((1 + z) / 3.2) ** 6.2)
    k = 1
    dV_dz = Distance_C(z)**2 / np.sqrt( omega_m*(1+z)**3 + 1 - omega_m )

    '''
    PL --- γ = -1.1 #(-2.5, 2.5)#
    CSFH --- zc = 2.8 #(0.1, 8.)#
    CPL --- γ = -0.6  #(-2.5, 2.5)#  zc = 5.5  #(0.1, 8.)#
    TSE --- γ1 = -0.7 #(-2., 2.)#  γ2 = 1.1 #(-3., 3.)#   s = 1.9 #(0.1, 4.)#
    TSRD --- a = 1.3 #(-2., 3.)#    b = 5.8 #(-3., 6.)#  c = 4.2 #(0.1, 5.)#
    '''

    if model == 'PL':
        gamma = args[0]
        k = (1+z)**gamma

    if model == 'CPL':
        gamma = args[0]
        zc = args[1]
        k = (1+z)**gamma * np.exp(-z/zc)

    if model == 'CSFH':
        zc = args[0]
        k = np.exp(-z/zc)

    if model == 'TSE':
        gamma1 = args[0]
        gamma2 = args[1]
        k_s = args[2]
        k = (1+z)**gamma1 / (1+((1+z)/k_s)**(gamma1+gamma2))

    if model == 'TSRD':
        k_a = args[0]
        k_b = args[1]
        k_c = args[2]
        k = (1 + z) ** k_a / (1 + ((1 + z) / k_c) ** (k_a + k_b)) / SFH

    return k*SFH*dV_dz/(1+z)

def pz_norm(z):
    norm = 1/fz.integrate(pz, z_min, z_max)
    ret = norm*pz(z)
    return ret



if __name__ == '__main__':

    plt.plot(z_space, pz(z_space, model='CPL', args = [-2.21, 4.51]))
    plt.show()
