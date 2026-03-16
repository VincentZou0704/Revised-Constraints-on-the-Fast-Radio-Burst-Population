import numpy as np
from FRBpopulation import FuncZou as fz
from FRBpopulation.setup import *
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.interpolate import interp2d
from scipy.optimize import minimize
import FRBpopulation.prepare.CHIME_Lin as CL
import FRBpopulation.MCMC.pzE as pzE
import FRBpopulation.prepare.Fluence as fce
from FRBpopulation.MCMC import Chi_analyze


model = 'PL'
theta = np.array([1.85, 42.64, 0.43, 5.64, -2.34])

######
bin = 30


def index_synx(num_all, x_axis_all):
    x_axis = []
    num = []
    for i in range(len(num_all)):
        if num_all[i] != 0:
            num.append(num_all[i])
            x_axis.append(x_axis_all[i])
    x_axis = np.array(x_axis)
    num = np.array(num)
    return num, x_axis


######  z plot
def zplot():
    plt.figure()

    num, boundary_z, _ = plt.hist(CL.z0, bins=bin, alpha=0.6, range=(0, z_max), color = '#CC99FF')
    x = boundary_z[:-1] + (boundary_z[1] - boundary_z[0])/2
    # num_cumulate = np.array([np.sum(num[i:bin]) for i in range(bin)])

    num_z = index_synx(num, x)[0]
    x_z = index_synx(num, x)[1]

    y_error_z = np.sqrt(num_z)

    plt.plot(z_space, pzE.pobs_z_norm(theta, z_space)*np.sum(np.diff(boundary_z)[0] * num_z), color = 'orange')
    plt.errorbar(x_z, num_z, yerr=y_error_z, fmt='o', markersize=3, color = 'red', alpha = 0.4)
    plt.ylim(ymin = 0)
    plt.xlabel('z with '+str(model)+' model')


######  logE plot
def Eplot():
    plt.figure()

    num, boundary_lgE, _ = plt.hist(CL.lgE0, bins=bin, alpha=0.6, range=(lgE_min, lgE_max), color = '#CC99FF')
    x = boundary_lgE[:-1] + (boundary_lgE[1] - boundary_lgE[0])/2

    num_lgE = index_synx(num, x)[0]
    x_lgE = index_synx(num, x)[1]

    y_error_lgE = np.sqrt(num_lgE)

    plt.plot(lgE_space, pzE.pobs_E_norm(theta, lgE_space)*np.sum(np.diff(boundary_lgE)[0] * num_lgE), color = 'orange')
    plt.errorbar(x_lgE, num_lgE, yerr=y_error_lgE, fmt='o', markersize=3, color = 'red', alpha = 0.4)
    plt.ylim(ymin = 0)
    plt.xlabel(r'logE with '+str(model)+' model')


###### logF_v plot
def Fplot():
    plt.figure()

    num, boundary_lgfv, _ = plt.hist(CL.lgfv0, bins=25, alpha=0.6, range=(lgfv_min, lgfv_max), color = '#CC99FF')
    x = boundary_lgfv[:-1] + (boundary_lgfv[1] - boundary_lgfv[0])/2

    num_lgfv = index_synx(num, x)[0]
    x_lgfv = index_synx(num, x)[1]

    y_error_lgfv = np.sqrt(num_lgfv)

    plt.plot(lgfv_space, fce.PDF_lgfv_obs(theta, lgfv_space)*np.sum(np.diff(boundary_lgfv)[0] * num_lgfv), color = 'orange')
    plt.errorbar(x_lgfv, num_lgfv, yerr=y_error_lgfv, fmt='o', markersize=3, color = 'red', alpha = 0.4)
    plt.ylim(ymin = 0)
    plt.xlabel(r'log$F_\nu$ with '+str(model)+' model')


if __name__ == '__main__':

    print(Chi_analyze.chi2_all(theta, CL.lgfv0, CL.lgE0, CL.z0))
    zplot()
    plt.title('parameters: ' + str(theta))

    Eplot()
    plt.title('parameters: ' + str(theta))

    Fplot()
    plt.title('parameters: '+str(theta))

    plt.show()

