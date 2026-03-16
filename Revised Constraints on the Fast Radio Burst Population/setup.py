from frb_mcmc import splinedata as sd
import numpy as np
import matplotlib.pyplot as plt


def Distance_C(z):      # comoving distance
    Dc = (c_photon/h0_P)*spldc(z)
    return Dc

def Distance_L(z):
    return Distance_C(z)*(1+z)

h0_P = 67.4
ob_P = 0.0493
f_IGM_p, alpha0, F0, sigma_host0, emu0= 0.83, 0., 0.2, 1, 100
a1 = 0.315
a2 = 1 - a1
omega_m = 0.315
c_photon = 2.99792458*10**8
gravity_c = 6.67408*1e-11
m_proton = 1.672621637*1e-27
coefficient = 1e-41*(21*c_photon/64/np.pi/gravity_c/m_proton/3.085677581467**2) # 1pc = 3.085677581467e16 m

ie = 1e-11

epsilon = 8.854187817 * 1e-12
h_planck = 6.62607015 * 1e-34
m_electron = 9.10956 * 1e-31
q_electron = 1.602176634 * 1e-19

spldc = sd.splinedcz(a1)
splz = sd.splinehez(a1)
splc = sd.splinec0()
spla = sd.splineA()

scale = 300
Vc = 400    # central frequency of the telescope
co1 = 1e22 / 3.085677581467e22 ** 2  # modified coefficient in Fluence Function
z_min, z_max = 0.01, 3.
z_space = np.linspace(z_min, z_max, scale)

lgE_min, lgE_max = 37., 43.
lgE_space = np.linspace(lgE_min, lgE_max, scale)

E_alpha0, lgEc0= 1.9, 41.5

lgfv_min, lgfv_max = -0.5, 2.
lgfv_space = np.linspace(lgfv_min, lgfv_max, scale)

lgfv_max_th0 = 0.75
n0 = 3

ks = True
model_z = 'PL'
type = 'Gold'
truths = [1.85, 42.64, 0.43, 5.64, -2.34]
BIN = False
bin = 30

if __name__ == '__main__':
    print(coefficient)