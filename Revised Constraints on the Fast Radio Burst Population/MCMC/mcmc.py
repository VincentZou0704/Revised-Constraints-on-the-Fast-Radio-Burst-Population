from FRBpopulation.setup import *
import FRBpopulation.prepare.CHIME_Lin as CL
import xlwt
from FRBpopulation.MCMC import Chi_analyze
import emcee
import corner


def prior(theta, modelx):

    E_alpha, lgEc, lgfv_max_th, n, *model_args = theta

    #-0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39 < lgEc < 44 and -10 < gamma < 10 and 0.1 < zc < 10:
    if modelx == 'SFH':
        if -0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39. < lgEc < 44.\
                :
            return 0
    elif modelx == 'CSFH':
        if -0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39. < lgEc < 44.\
                and 0.1 < model_args[0] < 10:
            return 0
    elif modelx == 'PL':
        if -0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39. < lgEc < 44.\
                and -10 < model_args[0] < 10:
            return 0
    elif modelx == 'CPL':
        if -0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39. < lgEc < 44.\
                and -10 < model_args[0] < 10 and 0.1 < model_args[1] < 10:
            return 0
    else:
        if -0.49 < lgfv_max_th < 2. and 1 < n < 8 and 1 < E_alpha < 5 and 39. < lgEc < 44.\
                and -10 < model_args[0] < 10 and -10 < model_args[1] < 10 and 0.1 < model_args[2] < 10:
            return 0
    return -np.inf


def log_prob(theta, lgfv, lgE, z):
    lp = prior(theta, modelx=model_z)
    if np.isfinite(lp):

        ######
        return lp - Chi_analyze.chi2_all(theta, lgfv, lgE, z, modelx = model_z)
    return -np.inf


# save mcmc results
def save_mcmc_result(flat_samples, path):
    my_workbook = xlwt.Workbook()
    sheet = my_workbook.add_sheet('mcmc_result')
    for i in range(flat_samples.shape[0]):
        for j in range(flat_samples.shape[1]):
            sheet.write(i, j, flat_samples[i][j])
    my_workbook.save(path)


def labels(modelx):
    if modelx == 'SFH':
        lb = [ r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$', r'n']
    elif modelx == 'CSFH':
        lb = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$', r'n', r'$z_c$']
    elif modelx == 'PL':
        lb = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$', r'n', r'$\gamma$']
    elif modelx == 'CPL':
        lb = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$', r'n', r'$\gamma$', '$z_c$']
    elif modelx == 'TSE':
        lb = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^max$', r'n', '$\gamma_1$', '$\gamma_2$', '$s$']
    elif modelx == 'TSRD':
        lb = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^max$', r'n', '$a$', '$b$', '$C$']
    else:
        lb = None
    return lb


if type == 'Full':
    lgfv0 = np.array(sorted(CL.lgfv_all))
    z0= np.array(sorted(CL.z_all))
    lgE0 = np.array(sorted(CL.lgE_all))
else:
    lgfv0 = np.array(sorted(CL.lgfv0))
    z0 = np.array(sorted(CL.z0))
    lgE0 = np.array(sorted(CL.lgE0))

if __name__ == '__main__':

    # run mcmc

    ndim = len(labels(model_z))
    nwalkers = 5*ndim

    #######

    pos = truths + 1e-4 * np.random.randn(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args = [lgfv0, lgE0, z0])
    sampler.run_mcmc(pos, 2000, progress=True)

    flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)

    #######
    save_mcmc_result(flat_samples, r'F:\pythonProject1\process\mcmc\_' + type + '_'+ model_z + '.xlsx')

    theta = []
    for i in range(flat_samples.shape[1]):
        theta.append(np.percentile(flat_samples[:, i], 50))
    theta = np.array(theta)
    ret_chi2 = Chi_analyze.chi2_all(theta, lgfv0, lgE0, z0, model_z)
    print(model_z + '_'+ type + 'parameters:', np.round(theta, 2), '\n', 'chi2:', np.round(ret_chi2, 2))

    plt.figure()
    fig = corner.corner(
        flat_samples, labels=labels(model_z),
        quantiles=[0.1587, 0.5, 0.8413], show_titles=True, title_kwargs={"fontsize": 12}, smooth = 1, smooth1d = 1
    )

    #######
    plt.savefig(r'F:\pythonProject1\process\mcmc\_' + type + '_'+ model_z + '.png')
    plt.show()