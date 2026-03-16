import matplotlib.pyplot as plt
from FRBpopulation.setup import *


mu_w = 5.49
sigma_w = 2.46
def w_inc(w, mu_w, sigma_w):
    x = -0.5*((np.log(w/mu_w))/np.log(sigma_w))**2
    return  np.exp(x) /(np.sqrt(2*np.pi)*np.log(sigma_w)*w)

def E_fluence(z, alpha, F, delta_v):
    left = 4*np.pi*(c_photon/spldc(z)/(h0_P * 3.085677581467e-13))**2/(1+z)**alpha
    return left* delta_v * F * 1e-29 # (with a unit 'erg', 1e-7)

# print(E_fluence(1,-1.5,1,1e9))


def Phi_z(Phi_0, z, n):
    def SFR_z (z):
        return 1.0025738*(1+z)**2.7/(1+((1+z)/2.9)**5.6)
    return (Phi_0/(1+z)) * (SFR_z(z)/SFR_z(0))**n

def DV_Domega_Dz(z):
    DH = c_photon/h0_P
    DC = DH * spldc(z)
    DA = DC/(1+z)
    Ez = np.sqrt(omega_m * (1 + z) ** 3 + 1 - omega_m)
    return DH*(1+z)**2*DA**2/Ez

z = np.linspace(1e-2,3,1000)
plt.plot(z,DV_Domega_Dz(z)*Phi_z(1,z,1.7))
plt.show()