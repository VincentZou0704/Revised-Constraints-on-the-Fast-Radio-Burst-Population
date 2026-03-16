import FRBpopulation.prepare.CHIME_Lin as CL
from FRBpopulation.MCMC import Chi_analyze
import numpy as np


i = 0  # 0 for full sample, 1 for gold sample

type = ['Full', 'Gold']

if i == 0:
    lgfv0 = np.array(sorted(CL.lgfv_all))
    z0= np.array(sorted(CL.z_all))
    lgE0 = np.array(sorted(CL.lgE_all))
else:
    lgfv0 = np.array(sorted(CL.lgfv0))
    z0 = np.array(sorted(CL.z0))
    lgE0 = np.array(sorted(CL.lgE0))



theta_SFH = np.array([[ 2.12, 42.38, 0.42, 4.29], [2.11, 42.49, 0.42, 6.63]])
theta_PL = np.array([[1.82, 42.49, 0.33, 3.90, -2.65], [1.85, 42.64, 0.43, 5.64,-2.34]])
theta_CSFH = np.array([[1.83, 42.35, 0.45, 2.86, 0.66], [1.90, 42.71, 0.43, 5.69, 0.82]])
theta_CPL = np.array([[1.89, 42.51, 0.39, 3.79, 1.10, 0.51], [1.93, 42.78, 0.44, 5.76, 1.85, 0.46]])
theta_TSE = np.array([[1.97, 42.48, 0.43, 3.79, 6.26, 2.61, 1.14], [1.96, 42.76, 0.44, 5.80, 1.48, 3.25, 1.39]])
theta_TSRD = np.array([[1.96, 42.48, 0.42, 3.83, 5.91, 1.02, 1.29], [1.96, 42.74, 0.44, 5.79, 3.30, 2.11, 1.66]])

kx = [4, 5, 5, 6, 7, 7]
model_dict = {'SFH': [theta_SFH, kx[0]], 'PL': [theta_PL, kx[1]],
              'CSFH':[theta_CSFH, kx[2]], 'CPL': [theta_CPL, kx[3]],
              'TSE': [theta_TSE, kx[4]], 'TSRD': [theta_TSRD, kx[5]]}

def likelihood(theta, lgfv, lgE, z, modelx):
    return np.exp(-1/2*Chi_analyze.chi2_all(theta, lgfv, lgE, z, modelx))

def aic(theta, lgfv, lgE, z, modelx, k):
    return -2*np.log(likelihood(theta, lgfv, lgE, z, modelx)) + 2*k

def bic(theta, lgfv, lgE, z, modelx, k):
    return -2*np.log(likelihood(theta, lgfv, lgE, z, modelx)) + k*np.log(85)

if __name__ == '__main__':
    print(type[i])
    for key,value in model_dict.items():
        ret1 = aic(value[0][i], lgfv0, lgE0, z0, key, value[1])
        ret2 = bic(value[0][i], lgfv0, lgE0, z0, key, value[1])
        print(key,'\t', 'AIC test:',  np.round(ret1,2), '\t','BIC tset', np.round(ret2,2))



