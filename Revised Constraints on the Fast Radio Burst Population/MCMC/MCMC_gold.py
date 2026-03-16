from FRBpopulation.MCMC.mcmc import *

if __name__ == '__main__':

    lgfv0, lgfv_all = np.array(sorted(CL.lgfv0)), np.array(sorted(CL.lgfv_all))
    z0, z_all = np.array(sorted(CL.z0)), np.array(sorted(CL.z_all))
    lgE0, lgE_all = np.array(sorted(CL.lgE0)), np.array(sorted(CL.lgE_all))

    # run mcmc

    ####### labels = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^max$', r'n', '$z_c$','$\gamma_1$', '$\gamma_2$', '$s$']
    labels = [r'$\alpha$', '$E_c$', r'log$F_{\nu,th}^{max}$', r'n', '$a$', '$b$', '$C$']

    ndim = len(labels)
    nwalkers = 30

    #######
    truths = [ 2.11, 42.36,  0.57,  3.59,  8.9,   0.26,  1.27]
    pos = truths + 1e-4 * np.random.randn(nwalkers, ndim)

    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=[lgfv_all, lgE_all, z_all])
    sampler.run_mcmc(pos, 2000, progress=True)

    flat_samples = sampler.get_chain(discard=100, thin=10, flat=True)

    #######
    save_mcmc_result(flat_samples, r'F:\pythonProject1\process\mcmc\TSRD_full.xlsx')

    theta = []
    for i in range(flat_samples.shape[1]):
        theta.append(np.percentile(flat_samples[:, i], 50))
    theta = np.array(theta)
    ret_chi2 = Chi_analyze.chi2_all(theta, lgfv0, lgE0, z0)
    print('TSRD full parameters:', np.round(theta, 2), '\n', 'chi2:', np.round(ret_chi2, 2))

    plt.figure()
    fig = corner.corner(
        flat_samples, labels=labels,
        quantiles=[0.1587, 0.5, 0.8413], show_titles=True, title_kwargs={"fontsize": 12}, smooth=1, smooth1d=1
    )

    #######
    plt.savefig(r"F:\pythonProject1\process\figures\TSRD_full.png")
    plt.show()