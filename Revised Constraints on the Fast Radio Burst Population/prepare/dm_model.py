import numpy as np
import matplotlib.pyplot as plt
from FRBpopulation.setup import *
from FRBpopulation import FuncZou
from scipy.interpolate import InterpolatedUnivariateSpline as IUS




def dm_cosmic_average(z, f, alpha):     #  f is parameterized as f = f_IGM,0 * (1 + alpha * z/(1+z) )
    #f = 0.84
    #alpha = 0
    return coefficient*ob_P*h0_P*f * splz(z)*(1+alpha*z/(1+z))



def p_host(dm_host, sigma_host, emu):       # log_normal distribution of dm_host
    x = np.log(np.abs(dm_host/emu)+ie)**2/(2*sigma_host**2)
    return  np.exp(-x)/(np.sqrt(2*np.pi)*dm_host*sigma_host)


def p_cosmic(dm_host, dm_frb, z, f, alpha, F):      # Macquart relation p(△)
    dm_cosmic = dm_frb - dm_host/(1+z)
    delta = (dm_cosmic/dm_cosmic_average(z, f , alpha) + ie)
    sigma = np.abs(F/np.sqrt(z))
    c = splc(sigma)
    A = spla(sigma)
    x = (delta**(-3)-c)**2/18/sigma**2
    return  A*delta**(-3)*np.exp(- x)


def intepart(dm_host, sigma_host, emu, dm_frb, z, f, alpha, F):
    return p_host(dm_host, sigma_host, emu)*p_cosmic(dm_host, dm_frb, z, f, alpha, F)


def p_dm_at_z(dm_frb, z, sigma_host, emu, f, alpha, F):
    int1 = FuncZou.integrate(intepart, ie, dm_frb * (1 + z) - ie, sigma_host, emu, dm_frb, z, f, alpha, F)
    return int1


def p_z_at_dm(z, dm_frb, sigma_host, emu, f, alpha, F):
    z_segment = np.linspace(0.001, 3, 1000)
    pz = FuncZou.func_build(z, z_segment, p_dm_at_z(dm_frb, z_segment, sigma_host, emu, f, alpha, F), mode = 'sline')
    return pz



#    # print(p_z_at_dm(np.linspace(0.001,2,1000), 1000, sigma_host0, emu0, f_IGM_p, alpha0, F0))
#    dm0 =1000
#    print(FuncZou.integrate(p_z_at_dm, 0.001, np.linspace(0.001,3,100), 1000, dm0, sigma_host0, emu0, f_IGM_p, alpha0, F0))