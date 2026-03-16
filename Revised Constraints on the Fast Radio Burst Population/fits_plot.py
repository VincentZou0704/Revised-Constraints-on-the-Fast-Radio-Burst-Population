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


i = 1  # 0 for full sample, 1 for gold sample

type = ['Full', 'Gold']

theta_SFH = np.array([[ 2.12, 42.38, 0.42, 4.29], [2.11, 42.49, 0.42, 6.63]])
theta_PL = np.array([[1.82, 42.49, 0.33, 3.90, -2.65], [1.85, 42.64, 0.43, 5.64,-2.34]])
theta_CSFH = np.array([[1.83, 42.35, 0.45, 2.86, 0.66], [1.90, 42.71, 0.43, 5.69, 0.82]])
theta_CPL = np.array([[1.89, 42.51, 0.39, 3.79, 1.36, 0.47], [1.94, 42.76, 0.44, 5.76, 2.17, 0.43]])
theta_TSE = np.array([[1.90, 42.38, 0.45, 3.08, 7.33, 2.86, 1.10], [1.92, 42.50, 0.43, 5.35, 0.95, 3.32, 1.42]])
theta_TSRD = np.array([[1.90, 42.34, 0.45, 3.09, 7.71, 0.95, 1.18], [1.92, 42.51, 0.43, 5.36, 2.99, 2.08, 1.69]])

colors = [ '#FF66FF', '#6699FF', '#66FF66', '#FF9933', '#FF3300', '#FFFF00']
model_dict = {'SFH': [theta_SFH, colors[0]], 'PL': [theta_PL, colors[1]],
              'CSFH':[theta_CSFH, colors[2]], 'CPL': [theta_CPL, colors[3]],
              'TSE': [theta_TSE, colors[4]], 'TSRD': [theta_TSRD, colors[5]]}

color_error = ['#6600FF', '#CC0000']

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
    z0 = [CL.z_all, CL.z0]
    chi2z = []

    num, boundary_z = np.histogram(z0[i], bins=bin, range=(0, z_max))
    x = boundary_z[:-1] + (boundary_z[1] - boundary_z[0])/2
    # num_cumulate = np.array([np.sum(num[i:bin]) for i in range(bin)])

    num_z = index_synx(num, x)[0]
    x_z = index_synx(num, x)[1]

    y_error_z = np.sqrt(num_z)
    for key, value in model_dict.items():
        chi2z.append(Chi_analyze.chi2_z(value[0][i], z0[i], key))
        plt.plot(z_space,
                 pzE.pobs_z_norm(value[0][i], z_space, modelx = key)*np.sum(np.diff(boundary_z)[0] * num_z),
                 color = value[1],
                 label = key)
    plt.errorbar(x_z, num_z, yerr=y_error_z, fmt='o', markersize=3,
                 color=color_error[i], alpha=0.7, label=type[i])

    plt.ylim(ymin = 0)
    plt.xlabel(r'$z$')
    plt.ylabel(r'N($z$) of ' + type[i] + ' sample')
    plt.legend(loc = 'upper right')
    return chi2z


######  logE plot
def Eplot():
    plt.figure()
    lgE0 = [CL.lgE_all, CL.lgE0]
    chi2lgE = []

    num, boundary_lgE = np.histogram(lgE0[i], bins=bin, range=(lgE_min, lgE_max))
    x = boundary_lgE[:-1] + (boundary_lgE[1] - boundary_lgE[0])/2

    num_lgE = index_synx(num, x)[0]
    x_lgE = index_synx(num, x)[1]

    y_error_lgE = np.sqrt(num_lgE)
    for key, value in model_dict.items():
        chi2lgE.append(Chi_analyze.chi2_lgE(value[0][i], lgE0[i], key))
        plt.plot(lgE_space,
                 pzE.pobs_E_norm(value[0][i], lgE_space, modelx = key)*np.sum(np.diff(boundary_lgE)[0] * num_lgE),
                 color=value[1],
                 label=key)
    plt.errorbar(x_lgE, num_lgE, yerr=y_error_lgE, fmt='o', markersize=3,
                 color=color_error[i], label=type[i], alpha=0.7)

    plt.ylim(ymin = 0)
    plt.xlabel(r'log$E$')
    plt.ylabel(r'N(log$E$) of ' + type[i] + ' sample')
    plt.legend(loc = 'upper right')

    return chi2lgE


###### logF_v plot
def Fplot():
    plt.figure()
    lgfv0 = [CL.lgfv_all, CL.lgfv0]
    chi2lgfv = []

    num, boundary_lgfv = np.histogram(lgfv0[i], bins=25, range=(lgfv_min, lgfv_max))
    x = boundary_lgfv[:-1] + (boundary_lgfv[1] - boundary_lgfv[0])/2

    num_lgfv = index_synx(num, x)[0]
    x_lgfv = index_synx(num, x)[1]

    y_error_lgfv = np.sqrt(num_lgfv)
    for key, value in model_dict.items():
        chi2lgfv.append(Chi_analyze.chi2_lgfv(value[0][i], lgfv0[i], key))

        #smooth
        window_size = 5
        # Calculate the moving average
        y = fce.PDF_lgfv(value[0][i], lgfv_space, modelx = key)*np.sum(np.diff(boundary_lgfv)[0] * num_lgfv)
        smoothed_y = np.convolve(y, np.ones(window_size) / window_size, mode='valid')

        plt.plot(lgfv_space[window_size-1:], smoothed_y,
                 color=value[1],
                 label=key)
    plt.errorbar(x_lgfv, num_lgfv, yerr=y_error_lgfv, fmt='o', markersize=3,
                 color=color_error[i], label=type[i], alpha = 0.7)

    plt.ylim(ymin = 0)
    plt.xlabel(r'log$F_\nu$')
    plt.ylabel(r'N(log$F_\nu$) of ' + type[i] + ' sample')
    plt.legend(loc = 'upper right')

    return chi2lgfv


chi2z = np.array(zplot())
# plt.savefig(r'F:\pythonProject1\figure\M1\fitsplot_' + type[i] + '_z.pdf')

chi2lgE = np.array(Eplot())
# plt.savefig(r'F:\pythonProject1\figure\M1\fitsplot_' + type[i] + '_E.pdf')

chi2lgfv = np.array(Fplot())
# plt.savefig(r'F:\pythonProject1\figure\M1\fitsplot_' + type[i] + '_Fv.pdf')

df = 85 - np.array([4, 5, 5, 6, 7, 7])
chi2all = (chi2z + chi2lgE + chi2lgfv) / df
print(np.round(chi2all, 2))

if __name__ == '__main__':
    plt.show()