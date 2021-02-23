# constants
import numpy as np
import math


N_D = 1e18 * 1e-12                       # Doper concentration [um^-3]
N_A = 0    * 1e-12                       # Acceptor concentration [um^-3]

e = 1.6e-19                             # charge [C]
eps0 = 8.85e-18                         # unit: C/(V*um)
m_0 = 9.1095e-31                         # unit: kg
T = 300                                 # Temperature [K]
kB = 1.38e-23                           # Boltzmann constant [m^2 kg s^-2 K^-1]
kBT_e = kB*T/e                          # unit: V
h = 6.62697e-34                         # unit:Js, Plank's constant
hv = 2.33                               # eV, energy per photon
sigma = 0.2                             # unit: um, SD of the gaussian distribution [um]

Az = 0.4                                # thickness [um^2]
V_SD = 0.00                             # unit: V
L_w = 10                                # unit: um, channel length of the nanowire
L_h = 10
Vbh_l = 0.3                             # unit: V, barrier height_left
Vbh_r = 0.3                             # unit: V, barrier height_right

# material
epsr = 30                               # unit: 1 'relative permittivity' (at 300K)   
eps_e = epsr*eps0/e                     # unit: 1/(V*um)
alpha = 1                               # 532 nm, absorption [um^-1]

t_n = 100 * 1e-9                        # unit: s          'lifetime of electron'
t_p = 100 * 1e-9                        # unit: s
t_s = 100 * 1e-12                       # unit: s
mu_n = 1000                             # unit: cm^2/(V*s)  'mobility of electron'
mu_p = 1000                             # unit: cm^2/(V*s), mobility of hole
m_e = 0.12 * m_0                        # unit: kg
m_h = 0.24 * m_0                        # unit: kg
Eg0_e = 0.3                             # unit: V,'band gap'
D_n = kBT_e * mu_n                      # diffusiion constant of e [cm^2/s]
D_p = kBT_e * mu_p                      # diffusiion constant of e [cm^2/s]

N_V = 2 * (2 * math.pi * m_h * kB * T / h**2)**(3/2)*1e-18      # unit: 1/um^3
N_C = 2 * (2 * math.pi * m_e * kB * T / h**2)**(3/2)*1e-18      # unit: 1/um^3
n_i = np.sqrt(N_C * N_V * np.exp(-Eg0_e/kBT_e))                 # unit: 1/um^3
EC_EF = kBT_e * np.log(N_C/N_D)

# Boundary Conditions
Psi_bi_l = Vbh_l - EC_EF                            # unit: V, build-in potential on left
Psi_bi_r = Vbh_r - EC_EF                            # unit: V, build-in potential on right
W_D_l = np.sqrt(2*eps_e/N_D*(Psi_bi_l - kBT_e))     # unit: um, depletion widtion
W_D_r = np.sqrt(2*eps_e/N_D*(Psi_bi_r - kBT_e))     # unit: um, depletion widtion
E_m_l = 2*(Psi_bi_l - kBT_e)/W_D_l                  # unit: V/um, maximum field at contact
E_m_r = 2*(Psi_bi_r - kBT_e)/W_D_r                  # unit: V/um, maximum field at contact

V_bc_l = -Vbh_l + V_SD                              # left BDC for V (bias applied on the left)
V_bc_r = -Vbh_r                                     # right BDC for V
n_up_bc_l = n_i*np.exp((Eg0_e/2 - Vbh_l)/kBT_e)     # unit: 1/um^3  
n_up_bc_r = n_i*np.exp((Eg0_e/2 - Vbh_r)/kBT_e)     # unit: 1/um^3
n_down_bc_l = n_i*np.exp((Eg0_e/2 - Vbh_l)/kBT_e)   # unit: 1/um^3  
n_down_bc_r = n_i*np.exp((Eg0_e/2 - Vbh_r)/kBT_e)   # unit: 1/um^3
p_bc_l = n_i*np.exp(-(Eg0_e/2 - Vbh_l)/kBT_e)       # unit: 1/um^3
p_bc_r = n_i*np.exp(-(Eg0_e/2 - Vbh_r)/kBT_e)       # unit: 1/um^3


# Scaled parameters
a = 10.0   # unit: um

N_D_bar = N_D/N_D
N_A_bar = N_A/N_D
n_i_bar = n_i/N_D

# sigma = 0.2um, SD of the gaussian distribution [um]
sigma_bar = sigma/a

# alpha = 1, absorption [um^-1]
alpha_bar = alpha * a

# hv = 2.33                             # eV, energy per photon
# mu_n = cm^2/Vs -> m^2/Vs, D_n: cm^2/s -> m^2/s

def D_factor(mu):
    return eps0 * epsr / ( mu * N_D * e * a**2 )
#     return (eps0*1e6) * epsr / ( (mu * 1e-4) * (N_D * 1e18) * e * (a*1e-6)**2 )
D_n_bar  = D_n * D_factor(mu_n)        # scaled value for D_n in the weak form of n
D_p_bar  = D_p * D_factor(mu_p)        # scaled value for D_p in the weak form of p

def t_p_factor(mu):
    # mu_n = 1000                             # unit: cm^2/(V*s)  'mobility of electron'
    # N_D = 1e18 * 1e-12                      # Doper concentration [um^-3]
    # eps0 = 8.85e-18                         # unit: C/(V*um)
    # t_n = 100 * 1e-9                        # unit: s          'lifetime of electron'
    # t_p = 100 * 1e-9                        # unit: s
    return mu * N_D * e / (eps0 * epsr) * 1e8
t_pn_bar = t_p * t_p_factor(mu_n)   # scaled value for t_p in the weak form of n
t_pp_bar = t_p * t_p_factor(mu_n)    # scaled value for t_p in the weak form of p

# G0 = 1/(um * s) -> 1/(1e-6 m *s), Az = um^2 -> m^2 * 1e-12, mu_n = cm^2/Vs -> 1e-4 * m^2/Vs
# a = um -> 1e-6 * m , eps0 = 8.85e-18                         # unit: C/(V*um)
V_scale_factor = (eps0 * epsr / (a**2 * e * N_D))
Vbh_l_bar = V_bc_l * V_scale_factor
Vbh_r_bar = V_bc_r * V_scale_factor
V_SD_bar = V_SD * V_scale_factor
