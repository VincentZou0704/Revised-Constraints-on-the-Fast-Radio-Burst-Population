from FRBpopulation.setup import *
from FRBpopulation import FuncZou
import pandas as pd
import emcee
import corner
import matplotlib.pyplot as plt


def modz(model, z, *args):
    SFH = (1+z)**2.6 / (1+((1+z)/3.2)**6.2)
    if model == '0':
        k = 1
    if model == 'PL':
        gamma = args[0]
        k = (1+z)**gamma
    if model == 'CPL':
        gamma = args[0]
        zc = args[1]
        k = (1+z)**gamma * np.exp(-z/zc)
    return k*SFH


def dm_cosmic_average(z, f, alpha):        #   f parameterized as f = f_IGM,0 * (1 + alpha * z/(1+z) )
    #f = 0.84
    #alpha = 0
    return coefficient*ob_P*h0_P*f * splz(z)*(1+alpha*z/(1+z))


def likelihood_host(dm_host, sigma_host, emu):
    x = np.log(np.abs(dm_host/emu)+ie)**2/(2*sigma_host**2)
    return  np.exp(-x)/(np.sqrt(2*np.pi)*dm_host*sigma_host)


def likelihood_cosmic(dm_host, dm_frb, z, f, alpha, F):
    dm_cosmic = dm_frb - dm_host/(1+z)
    delta = (dm_cosmic/dm_cosmic_average(z, f , alpha) + ie)
    sigma = np.abs(F/np.sqrt(z))       #fix F = 0.2
    c = splc(sigma)
    A = spla(sigma)
    x = (delta**(-3)-c)**2/18/sigma**2
    return  A*delta**(-3)*np.exp(- x)


def likelihood_DM(dm_host, sigma_host, emu, dm_frb, z, f, alpha, F):
    return likelihood_host(dm_host, sigma_host, emu)*likelihood_cosmic(dm_host, dm_frb, z, f, alpha, F)


def p_zDM(dm_frb, sigma_host, emu,  z, f, alpha, F, gamma, zc):
    int1 = FuncZou.integrate(likelihood_DM, ie, dm_frb*(1+z)-ie, 1000, sigma_host, emu, dm_frb, z, f, alpha, F)
    return modz('CPL', z, gamma, zc)*int1


def likelihood_zDM(theta, z, dm_frb):
    gamma, zc, F, sigma_host, emu = theta
    ret = np.array(p_zDM(dm_frb, sigma_host, emu, z, f_IGM_p, alpha0, F, gamma, zc)) + ie
    return np.sum(np.log(ret))


def log_prior(theta):
    gamma, zc, F, sigma_host, emu = theta
    if -2.5 < gamma < 2.5 and 0.1 < zc < 8 and 0.01 < F < 0.5 and 0.2 < sigma_host < 2 and 20 < emu < 200:
        return 0.0
    return -np.inf

def log_probability(theta, z, dm_frb):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + likelihood_zDM(theta, z, dm_frb)



def readsample(i):
    df = pd.read_excel(r'F:\pythonProject1\FRBpopulation\data\selected_samples.xlsx', usecols=[i])
    df_li = df.values.tolist()
    return np.array(df_li)

dm_frb, z0= readsample(0), readsample('z1')


truths = [-0.6, 5.5, F0, sigma_host0, emu0]
nwalkers, ndim =30, len(truths)
pos = truths + 1e-4*np.random.randn(nwalkers, ndim)
labels = [r'$\gamma$', r'$z_c$', 'F', '$\sigma_{host}$', '$e^\mu$']

sampler = emcee.EnsembleSampler(
    nwalkers, ndim, log_probability, args=(z0, dm_frb)
)

sampler.run_mcmc(pos, 1000, progress=True)

flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)
print(flat_samples.shape)

plt.figure()
fig = corner.corner(
    flat_samples, truths=truths ,labels=labels,
    quantiles=[0.1587, 0.5, 0.8413], show_titles=True, smooth = 1, smooth1d = 1
)
plt.savefig(r"F:\pythonProject1\FRBpopulation\progress\test.png")
plt.show()