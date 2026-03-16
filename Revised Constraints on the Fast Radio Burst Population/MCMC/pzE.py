import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
from scipy.optimize import minimize
from FRBpopulation.setup import *
import emcee
import FRBpopulation.prepare.CHIME_Lin as CL
import FRBpopulation.prepare.Fluence as fce
import xlwt
import corner
import FRBpopulation.prepare.z_model as zm
import FRBpopulation.prepare.selection_Fv as sf
from FRBpopulation import FuncZou as fz
from scipy.integrate import quad_vec
import time


###########
def pobs_lgE_z(lgE, z, theta, modelx):

    E_alpha, lgEc, lgfv_max_th, n, *model_args = theta

    lgfv = np.log10(fce.fluence(lgE, z))
    pz = zm.pz(z, model = modelx, args=model_args)
    f = pz*fce.p_lgE(lgE, E_alpha, lgEc)*sf.Ratio(lgfv, lgfv_max_th=lgfv_max_th, n=n)
    return f


def pobs_z(theta, z, modelx):
    integ = lambda lgE: pobs_lgE_z(lgE, z, theta, modelx)   # lgE_space and z_space are not the same shape

    # to calculate
    margin_E = fz.integrate_vec(integ, lgE_min, lgE_max)

    # to plot
    # margin_E = quad_vec(integ, lgE_min, lgE_max)[0]

    return margin_E


def pobs_z_norm(theta, z, modelx):
    integ_norm = lambda z: pobs_z(theta, z, modelx)
    norm = 1/fz.integrate(integ_norm, z_min, z_max)
    return norm*integ_norm(z)


def cdf_z_old(theta, z, modelx):
    integ = lambda z: pobs_z_norm(theta, z, modelx)
    ret = []
    for zi in z:
        ret.append(fz.integrate(integ, z_min, zi))
    return 1 - np.array(ret)


def cdf_z(theta, z, modelx):
    arr_a, arr_b = z_space, pobs_z(theta, z_space, modelx)
    # must give two array in advance
    integ = lambda z: fz.func_build(z, arr_a, arr_b)
    rets = fz.integrate(integ, z_min, z)
    ret =1 - rets/rets[-1]
    return ret



###########
def pobs_z_lgE(z, lgE, theta, modelx):

    E_alpha, lgEc, lgfv_max_th, n, *model_args= theta

    lgfv = np.log10(fce.fluence(lgE, z))
    pz = zm.pz(z, model = modelx, args=model_args)
    f = pz*fce.p_lgE(lgE, E_alpha, lgEc)**sf.Ratio(lgfv, lgfv_max_th=lgfv_max_th, n=n)
    return f


def pobs_E(theta, lgE, modelx):
    integ = lambda z: pobs_z_lgE(z, lgE, theta, modelx)
    # the definition of z_max doesn't matter
    # to calculate
    margin_z = fz.integrate_vec(integ, z_min, z_max)

    # to plot
    # margin_z = quad_vec(integ, z_min, z_max)[0]

    return margin_z


def pobs_E_norm(theta, lgE, modelx):
    integ_norm = lambda lgE: pobs_E(theta, lgE, modelx)
    norm = 1/fz.integrate(integ_norm, lgE_min, lgE_max)
    return norm*integ_norm(lgE)


def cdf_lgE_old(theta, lgE, modelx):
    integ = lambda lgE: pobs_E_norm(theta, lgE, modelx)
    ret = []
    for lgEi in lgE:
        ret.append(fz.integrate(integ, lgE_min, lgEi))
    return 1 - np.array(ret)


def cdf_lgE(theta, lgE, modelx):
    arr_a, arr_b = lgE_space, pobs_E(theta, lgE_space, modelx)
    # must give two array in advance
    integ = lambda lgE: fz.func_build(lgE, arr_a, arr_b)
    rets = fz.integrate(integ, lgE_min, lgE)
    ret =1 - rets/rets[-1]
    return ret


if __name__ == '__main__':
    lgfv0 = sorted(CL.lgfv0)
    z0 = sorted(CL.z0)
    lgE0 = sorted(CL.lgE0)
    theta0 = [1.9, 41.5]

    ###### test


    # y_ticks = [1e-1, 1, 10 , 1e2, 1e3]
    # y_axis = np.linspace(len(z0), 1, len(z0))
    # theta2 = np.array([ 2.30749667, 40.])
    # # [ 1.36489136, 40.] with p_int [ 2.30749667 40.] with p_obs
    # plt.plot(z_space, cdf_z(theta2, z_space)*(len(z0)))
    # plt.errorbar(z0, y_axis, yerr=np.sqrt(y_axis), fmt='o', markersize=3)
    # plt.yscale('log')
    # plt.ylim(ymin = y_ticks[0])
    # plt.xlabel('z')
    # plt.ylabel('N>z')
    # plt.savefig(r'F:\pythonProject1\FRBpopulation\progress\p_N_z.png')
    # plt.show()


    # y_ticks = [1e-1, 1, 10 , 1e2, 1e3]
    # y_axis = np.linspace(len(lgE0), 1, len(lgE0))
    # theta1 = np.array([ 2.07955802,43.])
    # # [ 1.65928797 42.14557502] with p_int [ 2.07955802,43.] with p_obs
    # plt.plot(lgE_space, cdf_lgE(theta1, lgE_space)*(len(lgE0)))
    # plt.errorbar(lgE0, y_axis, yerr=np.sqrt(y_axis), fmt='o', markersize=3)
    # plt.yscale('log')
    # plt.ylim(ymin = y_ticks[0])
    # plt.xlabel('logE')
    # plt.ylabel('N>lgE')
    # plt.savefig(r'F:\pythonProject1\FRBpopulation\progress\p_N_lgE.png')
    # plt.show()


    # plt.figure()
    # theta1 = np.array([1.8,41.5,-2,5])
    # plt.plot(lgE_space, pobs_E_norm(theta1, lgE_space), 'r')
    # plt.xlabel('logE')
    # plt.ylabel('p_value')
    # plt.title('p(logE) marginalized at z')
    # plt.figure()
    # theta2 = theta0
    # plt.plot(z_space, pobs_z_norm(theta1, z_space), 'g')
    # plt.xlabel('z')
    # plt.ylabel('p_value')
    # plt.title('p(z) marginalized at logE')
    # plt.show()

    truths = [1.85, 42.64, 0.43, 5.64, -2.34]
    start = time.time()
    print(cdf_lgE(truths,np.array(lgE0),'PL'))
    plt.figure()
    plt.plot(lgE_space, cdf_lgE(truths, lgE_space, 'PL'), 'r')
    plt.show()

    end = time.time()
    print('Running time: ', end - start, 's')





    ###### run


    # labels = [r'$\alpha$', '$E_c$']
    #
    # ndim = len(labels)
    # nwalkers = 5 * ndim
    # truths = [E_alpha0, lgEc0]
    # pos = truths + 1e-4 * np.random.randn(nwalkers, ndim)
    #
    # # [] or () should be used in args
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_E, args = [lgE0])
    # sampler.run_mcmc(pos, 1000, progress=True)
    #
    # flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)
    # save_mcmc_result(flat_samples)
    #
    # plt.figure()
    # fig = corner.corner(
    #     flat_samples, truths=truths ,labels=labels,
    #     quantiles=[0.1587, 0.5, 0.8413], show_titles=True, title_kwargs={"fontsize": 12}, smooth = 1, smooth1d = 1
    # )
    # plt.savefig(r"F:\pythonProject1\process\figures\pE_CHIME.png")
    # plt.show()

