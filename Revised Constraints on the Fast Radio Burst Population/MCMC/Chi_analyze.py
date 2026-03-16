from FRBpopulation.setup import *
import FRBpopulation.prepare.Fluence as fce
import FRBpopulation.prepare.CHIME_Lin as CL
from FRBpopulation.MCMC import pzE
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy.stats import chi2
import time




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


def chi2_lgE(theta, lgE, modelx):

    if BIN:
        num_all, boundary = np.histogram(lgE, bins=bin, range=(lgE_min, lgE_max))
        x_axis_all = boundary[:-1] + (boundary[1] - boundary[0])/2

        num = index_synx(num_all, x_axis_all)[0]
        x_axis = index_synx(num_all, x_axis_all)[1]

        y = pzE.pobs_E_norm(theta, x_axis, modelx)*np.sum(np.diff(boundary)[0] * num)
        y_axis = num
        sigma2 = num

    else:
        lgE = np.array(sorted(lgE))
        y_axis = np.linspace(len(lgE), 1., len(lgE))
        sigma2 = y_axis
        y = pzE.cdf_lgE(theta, lgE, modelx)*len(lgE)

    if ks:
        return y, y_axis

    chi2 = (y_axis - y)**2/sigma2

    return np.sum(chi2)


def chi2_z(theta, z, modelx):

    if BIN:
        num_all, boundary = np.histogram(z, bins=bin, range=(0, z_max))
        x_axis_all = boundary[:-1] + (boundary[1] - boundary[0]) / 2

        num = index_synx(num_all, x_axis_all)[0]
        x_axis = index_synx(num_all, x_axis_all)[1]

        y = pzE.pobs_z_norm(theta, x_axis, modelx) * np.sum(np.diff(boundary)[0] * num)
        y_axis = num
        sigma2 = num

    else:
        z = np.array(sorted(z))
        y_axis = np.linspace(len(z), 1., len(z))
        sigma2 = y_axis
        y = pzE.cdf_z(theta, z, modelx)*len(z)

    if ks:
        return y, y_axis
    chi2 = (y_axis - y)**2/sigma2
    ret = np.sum(chi2)
    # print(theta, ret)
    return ret


def chi2_lgfv(theta, lgfv, modelx):

    if BIN:

        num_all, boundary = np.histogram(lgfv, bins=25, range=(lgfv_min, lgfv_max))
        x_axis_all = boundary[:-1] + (boundary[1] - boundary[0]) / 2

        num = index_synx(num_all, x_axis_all)[0]
        x_axis = index_synx(num_all, x_axis_all)[1]

        ######
        y = fce.PDF_lgfv(theta, x_axis, modelx) * np.sum(np.diff(boundary)[0] * num)
        y_axis = num
        sigma2 = num

    else:
        lgfv = np.array(sorted(lgfv))
        y_axis = np.linspace(len(lgfv) , 1., len(lgfv))
        sigma2 = y_axis
        y = fce.CDF_lgfv_sf(theta, lgfv, modelx)*(len(lgfv))

    if ks:
        return y, y_axis

    chi2 = (y_axis - y)**2/sigma2
    ret = np.sum(chi2)
    # print(theta, ret)
    return ret


def chi2_all(theta, lgfv, lgE, z, modelx):
    return chi2_lgE(theta, lgE, modelx) + chi2_z(theta, z, modelx) + chi2_lgfv(theta, lgfv, modelx)


if __name__ == '__main__':



    start = time.time()
    lgfv0 = np.array(sorted(CL.lgfv_all))
    z0 = np.array(sorted(CL.z_all))
    lgE0 = np.array(sorted(CL.lgE_all))

    # (1, 5), (39, 44), (-0.49, 2), (1,8), (-10, 10), (0.1, 10)
    bounds = [(1, 5), (39, 44), (-0.49, 2), (1, 8),  (-10, 10), (-10, 10), (0.1, 10)]

    # 1.9, 41.5, 0.75, 3, -0.6, 5.5
    initial = np.array([1.9, 41.5, 0.75, 3, 1.3, 5.8, 4.2])


    ######Bin:


    # df = len(z0)
    # cl95_lower, cl68_lower, medium, cl68_upper, cl95_upper = chi2.ppf([0.025, 0.16, 0.5, 0.84, 0.975], df = df)
    # cl_chi2 = [cl95_lower, cl68_lower, medium, cl68_upper, cl95_upper]
    # print(cl_chi2)
    # print(chi2_lgfv([1.01,40.83,0.43], lgfv0))

    if BIN:
        df = 3*bin - len(initial)
    else:
        df = len(z0) - len(initial)


    # soln1 = minimize(chi2_lgE, initial, bounds = bounds, args=(lgE0))
    # print('chi2_lgE:', np.round(chi2_lgE(soln1.x, lgE0), 2),
    #       'chi/df:', np.round(chi2_lgE(soln1.x, lgE0) / df, 2),
    #       'parameters: ', np.round(soln1.x, 2), '\n')
    #
    # soln2 = minimize(chi2_z, initial, bounds = bounds, args=(z0))
    # print('chi2_z:', np.round(chi2_z(soln2.x, z0), 2),
    #       'chi/df:', np.round(chi2_z(soln2.x, z0) / df, 2),
    #       'parameters: ', np.round(soln2.x, 2), '\n')

    # soln3 = minimize(chi2_lgfv, initial, bounds = bounds, args=(lgfv0))
    # print('chi2_lgfv:',np.round(chi2_lgfv(soln3.x, lgfv0), 2),
    #       'chi/df:', np.round(chi2_lgfv(soln3.x, lgfv0) / df, 2),
    #       'parameters: ', np.round(soln3.x, 2), '\n')

    soln = minimize(chi2_all, initial, bounds = bounds, args=(lgfv0, lgE0, z0))
    print('chi2_all:',np.round(chi2_all(soln.x, lgfv0, lgE0, z0), 2),
          'chi/df:', np.round(chi2_all(soln.x, lgfv0, lgE0, z0)/df, 2),
          'parameters: ', np.round(soln.x, 2), '\n')

    end = time.time()
    print('running time:',end - start, 's')