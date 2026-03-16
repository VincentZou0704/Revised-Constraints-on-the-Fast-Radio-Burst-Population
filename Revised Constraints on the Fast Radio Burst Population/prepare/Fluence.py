from FRBpopulation.setup import *
from FRBpopulation import FuncZou as fz
import matplotlib.pyplot as plt
from FRBpopulation.prepare import z_model as zm
from scipy.optimize import minimize
from scipy.optimize import fsolve
from scipy import integrate
import FRBpopulation.prepare.selection_Fv as sf
import time
from scipy.interpolate import interp2d


def fluence(lgE, z):      # log_Fluence at energy and redshift
    DL2 = Distance_L(z)**2
    fv = co1 * np.sqrt(1+z)*10**lgE / (4*np.pi*DL2*Vc)
    return fv

# print(fluence(40,1))

# intrinsic energy distribution at d(lgE)
def p_lgE_norm(lgE, E_alpha, lgEc):
    func = lambda lgE : (10**lgE / 10**lgEc) ** (-E_alpha + 1) * np.exp(-10**lgE / 10**lgEc)
    norm = 1/fz.integrate(func, lgE_min, lgE_max)
    x = 10**lgE/10**lgEc
    return norm*x**(-E_alpha + 1) * np.exp(-x)


def p_lgE(log_E, E_alpha, log_Ec):
    x = 10 ** log_E / 10 ** log_Ec
    return x ** (-E_alpha + 1) * np.exp(-x)


def iso_lgE(lgfv, z):    # isotropic energy at log_fluence and redshift
    DL2 = Distance_L(z)**2
    fe = (4 * np.pi * DL2 * Vc) / co1 * np.sqrt(1+z) * 10**lgfv
    return np.log10(fe)


def CDF_lgfv_sf(theta, lgfv, modelx):
    # CDF of fluence with selection function (with parameters)

    E_alpha, lgEc, lgfv_max_th, n, *model_args = theta

    lgE_upperlimit = lambda z: np.log10(10**lgfv/(co1 * np.sqrt(1+z) / (4*np.pi*Distance_L(z)**2*Vc)))

    ######
    select = lambda lgE, z: sf.Ratio(np.log10(fluence(lgE, z)), lgfv_max_th = lgfv_max_th, n = n)

    pE = lambda lgE: p_lgE(lgE, E_alpha, lgEc)

    pz = lambda z:  zm.pz(z, model = modelx, args = model_args)
    ######

    integ_1 = lambda lgE, z:  select(lgE, z)* pE(lgE)

    integ_2 = lambda z :  pz(z) * fz.integrate(integ_1, lgE_min, lgE_upperlimit(z), z)

    up_max = lambda z: np.log10(10**lgfv_max/(co1 * np.sqrt(1+z) / (4*np.pi*Distance_L(z)**2*Vc)))
    up_min = lambda z: np.log10(10**lgfv_min/(co1 * np.sqrt(1+z) / (4*np.pi*Distance_L(z)**2*Vc)))

    max_integ = lambda z :  pz(z) * fz.integrate(integ_1, lgE_min, up_max(z), z)
    min_integ = lambda z :  pz(z) * fz.integrate(integ_1, lgE_min, up_min(z), z)

    cdf_max = fz.integrate(max_integ, z_min, z_max)
    cdf_min = fz.integrate(min_integ, z_min, z_max)

    # quad can save computation power while fz.fanc speed up by 1.5 time
    rets = fz.integrate_vec(integ_2, z_min, z_max)
    ret = np.clip(1 - (rets - cdf_min)/(cdf_max - cdf_min), 0, 1)

    # This term should be adapted to the usage
    return ret



def PDF_lgfv(theta, lgfv, modelx):
    dx = 1e-3
    dy = (CDF_lgfv_sf(theta, (lgfv - dx), modelx)- CDF_lgfv_sf(theta, lgfv, modelx))
    pdf = dy/dx

    # pfv = np.gradient( 1 - CDF_lgfv_sf(theta, lgfv_space), lgfv_space)
    # pdf = fz.func_build(lgfv, lgfv_space, pfv)
    return pdf




if __name__ == '__main__':

    import FRBpopulation.prepare.CHIME_Lin as CL

    theta0 = np.array([1.8, 41.5])
    fv0 = 10**np.array(sorted(CL.lgfv0))

    # soln = minimize(chi2_fv, theta0, bounds=([0,3],[40,43]), args=fv0)
    # print(soln.x)


    # theta3 = np.array([ 1.02764158, 40.])
    # # [ 1.02764158 40.] with p_int [2.36572553, 43.] with p_obs
    # y_ticks = [1e-1, 1, 10 , 1e2, 1e3]
    # y_axis = np.linspace(len(fv0), 1, len(fv0))
    # plt.plot(lgfv_space, CDF_fv(theta3, 10**lgfv_space)*(len(fv0)),'r')
    # plt.errorbar(np.log10(fv0), y_axis, yerr=np.sqrt(y_axis), fmt='.b', markersize=3)
    # plt.yscale('log')
    # plt.ylim(ymin = y_ticks[0])
    # plt.xlabel(r'log$F_\nu$')
    # plt.ylabel(r'N>log$F_\nu$')
    # # plt.savefig(r'F:\pythonProject1\FRBpopulation\progress\p_N_lgfv.png')
    # plt.show()


    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 1.0, 41.5), label=r'$\alpha = 1.0$', color='blue')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 1.2, 41.5), label=r'$\alpha = 1.2$', color='yellow')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 1.5, 41.5), label=r'$\alpha = 1.5$', color='black')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 1.8, 41.5), label=r'$\alpha = 1.8$', color='red')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 2.0, 41.5), label=r'$\alpha = 2.0$', color='orange')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 2.2, 41.5), label=r'$\alpha = 2.2$', color='green')
    # plt.plot(lgE_space, p_lgE_norm(lgE_space, 3.0, 41.5), label=r'$\alpha = 3.0$', color='purple')
    # plt.xlabel(r'logE, log$E_c$=41.5')
    # plt.ylabel(r'dN/dlogE')
    # plt.legend()
    # plt.show()


    start_time = time.time()

    # print(CDF_lgfv_sf((E_alpha0, lgEc0, 0.75, 3), np.linspace(-0.5,2,50)))

    theta = np.array([2.12, 42.64, 0.43, 5.64])
    plt.figure()
    # lgfv_space = np.linspace(0.3, 0.5, 100)
    plt.plot(lgfv_space, PDF_lgfv(theta, lgfv_space, modelx='SFH'), color = '#22DE4F', label = r'log$F_{\nu,obs}$')
    plt.xlabel(r'$\mathrm{log}F_\nu$')
    plt.ylabel('p')
    plt.title('Fluence distribution')
    plt.legend(loc='upper right')
    plt.show()

    end_time = time.time()
    print('Running time: ', end_time - start_time, 's')
