import sys
sys.path.insert(0, "/home/quantumer/dmde_class_perturb-CEDE/python/build/lib.linux-x86_64-3.10")
import subprocess

import classy
print("classy imported from: \n",classy.__file__)

# import necessary modules
# uncomment to get plots displayed in notebook
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import os
import math
from scipy.interpolate import UnivariateSpline

basedir = os.path.abspath('./')
gdmdir = os.path.join(basedir,'gdm_files')

gdm_exists = os.path.exists('./gdm_files')

if gdm_exists :
     print('gdm exists')
#    subprocess.run(["rm","-r","gdm_files"])
#    subprocess.run(["mkdir","gdm_files"])
else :
    subprocess.run(["mkdir","gdm_files"])

#############################################
#
# Cosmological parameters:
z_max_pk = 5e7       # highest redshift involved
z_max_rec = 1e6
k_per_decade =250     # number of k values, controls final resolution
k_min_tau0 = 0.1       # this value controls the minimum k value in the figure (it is k_min * tau0)
P_k_max_inv_Mpc =5   # this value is directly the maximum k value in the figure in Mpc
tau_num_early = 4200   # number of conformal time values before recombination, controls final resolution
tau_num_late = 250     # number of conformal time values after recombination, controls final resolution
tau_ini = 0.01          # first value of conformal time in Mpc
#tau_label_Hubble = 20. # value of time at which we want to place the label on Hubble crossing
#tau_label_ks = 40.     # value of time at which we want to place the label on sound horizon crossing
#tau_label_kd = 230.    # value of time at which we want to place the label on damping scale crossing
h = 0.6732117
#
# Cosmological parameters and other CLASS parameters
#
common_settings = {# which output? transfer functions only
                   'output':'mTk',
                   # LambdaCDM parameters
                   'h':0.6732117,
                   'omega_b':0.02238280,
                   'omega_cdm':0.12038,
                   'frac_dmde': 1.0,
                   'g_scf': -0.013,
                   'm_scf': 14.3,
                   'phi_scf_ini': 0.695,
                   'omega_ncdm': 0.0006451439,
                   'N_ncdm': 1,
                   'm_ncdm': 0.06,
                   'T_ncdm': 0.71611,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'Omega_k': 0,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.2454006,
                   # other output and precision parameters
                   'z_max_pk':z_max_pk,
                   # to get a larger z_max in recfast, 
                   # we must increase both recfast_z_initial 
                   # and the number of sampled values recfast_Nz0
                   # (in order to keep the same stepzize asd in the default: Delta z = 0.5)
                   'recfast_Nz0':z_max_rec*2.,
                   'recfast_z_initial':z_max_rec,
                   #'k_step_sub':'0.01',
                   'k_per_decade_for_pk':k_per_decade,
                   'k_per_decade_for_bao':k_per_decade,
                   'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                   'perturb_sampling_stepsize':'0.005',
                   'P_k_max_1/Mpc':P_k_max_inv_Mpc,
                   #'compute damping scale':'yes', # needed to output and plot Silk damping scale
                   'gauge':'synchronous'}

###############
#
# call CLASS
#
###############
print("Starting CLASS\n")

M_pre = Class()
M_pre.set(common_settings)
M_pre.compute()

print("CLASS is done!\n")

background = M_pre.get_background() # load background table                                                                                                                           
#print background.viewkeys()                                                                                                                                                          
background_tau = background['conf. time [Mpc]'] # read conformal times in background table       
background_z = background['z'] # read redshift
background_tau_at_z = interp1d(background_z,background_tau)

tau_ini = background_tau_at_z(z_max_pk)
print(tau_ini)
times = M_pre.get_current_derived_parameters(['tau_rec','conformal_age'])
tau_rec=times['tau_rec']
tau_0 = times['conformal_age']
tau1 = np.logspace(math.log10(tau_ini),math.log10(5*tau_rec),tau_num_early)
tau2 = np.logspace(math.log10(5*tau_rec),math.log10(tau_0),tau_num_late)[1:]
tau2[-1]=tau_0
tau2[-1] *= 0.999999 # this tiny shift avoids interpolation errors
tau= np.concatenate((tau1,tau2))
tau_num = len(tau)

background_z = background['z'] # read redshift
background_aH = background['H [1/Mpc]']/(1.+background['z'])

background_rho_cdm = background['(.)rho_cdm'] # read rho cdm in background table
background_rho_nu = background['(.)rho_ur'] # read rho nu in background table

background_z_at_tau = interp1d(background_tau,background_z)
a_tau = 1/(1+background_z_at_tau(tau))

background_aH_at_a = interp1d(1/(background_z+1),background_aH)
background_aH_at_tau = interp1d(background_tau,background_aH)
background_a_at_tau = interp1d(background_tau,1/(background_z+1))

background_rho_cdm_tau = interp1d(background_tau,background_rho_cdm)
background_rho_nu_tau = interp1d(background_tau,background_rho_nu)

background_rho_phi = background['(.)rho_scf']
background_p_phi = background['(.)p_scf']
background_rho_phi_tau = interp1d(background_tau,background_rho_phi)
background_p_phi_tau = interp1d(background_tau,background_p_phi)

background_rho_ncdm = background['(.)rho_ncdm[0]']
background_rho_ncdm_tau = interp1d(background_tau,background_rho_ncdm)
background_p_ncdm = background['(.)p_ncdm[0]']
background_p_ncdm_tau = interp1d(background_tau,background_p_ncdm)

# w_T_ncdm = background_p_ncdm/background_rho_ncdm
# w_T_phi = background_p_phi/background_rho_phi
w_T_general = (1/3*background_rho_nu+background_p_phi+background_p_ncdm)/(background_rho_cdm+background_rho_nu+background_rho_phi+background_rho_ncdm)
# w_T_general = (1/3*background_rho_nu+background_p_phi)/(background_rho_cdm+background_rho_nu+background_rho_phi)

w_T_at_a_interp = interp1d(1/(1+background_z),w_T_general)
w_T_at_a_list = w_T_at_a_interp(a_tau)

w_T_at_a_prime_interp = UnivariateSpline(1/(1+background_z),w_T_general,s=0).derivative()
w_T_at_a_prime_list = w_T_at_a_prime_interp(a_tau)

w_ncdm_at_a_interp = UnivariateSpline(1/(1+background_z),background_p_ncdm/background_rho_ncdm,s=0)
w_ncdm_at_a_list = w_ncdm_at_a_interp(a_tau)
w_ncdm_at_a_prime_interp = w_ncdm_at_a_interp.derivative()
w_ncdm_at_a_prime_list = w_ncdm_at_a_prime_interp(a_tau)
p_ncdm_at_a_interp = UnivariateSpline(1/(1+background_z),background_p_ncdm,s=0)
p_ncdm_at_a_list = p_ncdm_at_a_interp(a_tau)

p_phi_at_a_interp = UnivariateSpline(1/(1+background_z),background_p_phi,s=0)
p_phi_at_a_list = p_phi_at_a_interp(a_tau)
rho_phi_at_a_interp = UnivariateSpline(1/(1+background_z),background_rho_phi,s=0)
rho_phi_at_a_list = rho_phi_at_a_interp(a_tau)
w_phi_at_a_interp = UnivariateSpline(1/(1+background_z),background_p_phi/background_rho_phi,s=0)
w_phi_at_a_list = w_phi_at_a_interp(a_tau)


ones = np.ones(len(background_z))   
ca2_list = w_T_at_a_interp(1/(1+background_z)) - 1./3.*1/(1+background_z)*w_T_at_a_prime_interp(1/(1+background_z))/(ones+w_T_at_a_interp(1/(1+background_z)))
ca2_a = interp1d(1/(1+background_z),ca2_list)

ones = np.ones(len(a_tau))
ca2_list = w_T_at_a_list - 1./3.*a_tau*w_T_at_a_prime_list/(1+w_T_at_a_list)
ca2_tau = interp1d(tau,ca2_list)

w_T_tau = interp1d(background_tau,w_T_general)

print(background_z_at_tau(tau[0]))
print("Computing the perturbations\n")
one_time = M_pre.get_transfer(background_z_at_tau(tau[0]))
k = one_time['k (h/Mpc)']
k_num = len(k)

cs2_d = np.zeros((tau_num,k_num))
delta_d = np.zeros((tau_num,k_num))
theta_d = np.zeros((tau_num,k_num))
c2eff_d = np.zeros((tau_num,k_num))
Dc2_d = np.zeros((tau_num,k_num))

for i in range(tau_num):
    one_time = M_pre.get_transfer(background_z_at_tau(tau[i]))
    
    delta_d[i,:] = (one_time['d_ur'][:]*background_rho_nu_tau(tau[i]) + one_time['d_cdm'][:]*background_rho_cdm_tau(tau[i]) + one_time['d_scf'][:]*background_rho_phi_tau(tau[i]) + one_time['d_ncdm[0]'][:]*background_rho_ncdm_tau(tau[i]))\
    /(background_rho_cdm_tau(tau[i])+background_rho_nu_tau(tau[i])+background_rho_phi_tau(tau[i])+background_rho_ncdm_tau(tau[i]))
    
    cs2_d[i,:] = (1./3.*one_time['d_ur'][:]*background_rho_nu_tau(tau[i])+one_time['d_scf'][:]*background_p_phi_tau(tau[i])+one_time['d_ncdm[0]'][:]*background_p_ncdm_tau(tau[i]))\
            /(background_rho_cdm_tau(tau[i])*one_time['d_cdm'][:]+background_rho_nu_tau(tau[i])*one_time['d_ur'][:]+background_rho_phi_tau(tau[i])*one_time['d_scf'][:]+background_rho_ncdm_tau(tau[i])*one_time['d_ncdm[0]'][:])
    
    # delta_d[i,:] = (one_time['d_ur'][:]*background_rho_nu_tau(tau[i]) + one_time['d_cdm'][:]*background_rho_cdm_tau(tau[i]) + one_time['d_scf'][:]*background_rho_phi_tau(tau[i]))\
    # /(background_rho_cdm_tau(tau[i])+background_rho_nu_tau(tau[i])+background_rho_phi_tau(tau[i]))
    
    # cs2_d[i,:] = (1./3.*one_time['d_ur'][:]*background_rho_nu_tau(tau[i])+one_time['d_scf'][:]*background_p_phi_tau(tau[i]))\
    #         /(background_rho_cdm_tau(tau[i])*one_time['d_cdm'][:]+background_rho_nu_tau(tau[i])*one_time['d_ur'][:]+background_rho_phi_tau(tau[i])*one_time['d_scf'][:])
    
    # cs2_d[i,:] = (1./3.*one_time['d_ur'][:]*background_rho_nu_tau(tau[i]))\
    #          /(background_rho_cdm_tau(tau[i])*one_time['d_cdm'][:]+background_rho_nu_tau(tau[i])*one_time['d_ur'][:])
    
common_settings = {# which output? transfer functions only
                   'output':'vTk',
                   # LambdaCDM parameters
                   'h':0.6732117,
                   'omega_b':0.02238280,
                   'omega_cdm':0.12038,
                   'frac_dmde': 1.0,
                   'g_scf': -0.013,
                   'm_scf': 14.3,
                   'phi_scf_ini': 0.695,
                   'omega_ncdm': 0.0006451439,
                   'N_ncdm': 1,
                   'm_ncdm': 0.06,
                   'T_ncdm': 0.71611,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'Omega_k': 0,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.2454006,
                   # other output and precision parameters
                   'z_max_pk':z_max_pk,
                   # to get a larger z_max in recfast, 
                   # we must increase both recfast_z_initial 
                   # and the number of sampled values recfast_Nz0
                   # (in order to keep the same stepzize asd in the default: Delta z = 0.5)
                   'recfast_Nz0':z_max_rec*2.,
                   'recfast_z_initial':z_max_rec,
                   #'k_step_sub':'0.01',
                   'k_per_decade_for_pk':k_per_decade,
                   'k_per_decade_for_bao':k_per_decade,
                   'k_min_tau0':k_min_tau0, # this value controls the minimum k value in the figure
                   'perturb_sampling_stepsize':'0.005',
                   'P_k_max_1/Mpc':P_k_max_inv_Mpc,
                   #'compute damping scale':'yes', # needed to output and plot Silk damping scale
                   'gauge':'synchronous'}

M_pre = Class()
M_pre.set(common_settings)
M_pre.compute()

one_time = M_pre.get_transfer(background_z_at_tau(tau[0]))
k = one_time['k (h/Mpc)']
k_num = len(k)

for i in range(tau_num):
    one_time = M_pre.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau

#     delta_T_std[i]= (background_rho_cdm_at_a_std(a_std[i])*delta_cdm[i]+background_rho_nu_at_a_std(a_std[i])*delta_nu[i])\
#     /(background_rho_cdm_at_a_std(a_std[i])+background_rho_nu_at_a_std(a_std[i]))
#     theta_T_std[i] = ((1+1/3)*background_rho_nu_at_a_std(a_std[i])*theta_nu[i])\
#                     /((background_rho_cdm_at_a_std(a_std[i])+background_rho_nu_at_a_std(a_std[i])*(1.+1./3.)))

    theta_d[i,:] = ((1+1/3)*one_time['t_ur'][:]*background_rho_nu_tau(tau[i])+(1+w_phi_at_a_list[i])*one_time['t__scf'][:]*background_rho_phi_tau(tau[i])+(1+w_ncdm_at_a_list[i])*one_time['t_ncdm[0]']*background_rho_ncdm_tau(tau[i]))\
    /(background_rho_cdm_tau(tau[i])+(1+1/3)*background_rho_nu_tau(tau[i])+(1+w_phi_at_a_list[i])*background_rho_phi_tau(tau[i])+(1+w_ncdm_at_a_list[i])*background_rho_ncdm_tau(tau[i]))

    # theta_d[i,:] = ((1+1/3)*one_time['t_ur'][:]*background_rho_nu_tau(tau[i])+(1+w_phi_at_a_list(i))*one_time['t__scf'][:]*background_rho_phi_tau(tau[i]))\
    # /(background_rho_cdm_tau(tau[i])+(1+1/3)*background_rho_nu_tau(tau[i])+(1+w_phi_at_a_list(i))*background_rho_phi_tau(tau[i]))
    
    c2eff_d[i,:] = ((k*h)**2*cs2_d[i,:]*delta_d[i,:]+3*background_aH_at_tau(tau[i])*ca2_tau(tau[i])*(1+w_T_tau(tau[i]))*theta_d[i,:])\
        /((k*h)**2*delta_d[i,:]+3*background_aH_at_tau(tau[i])*(1+w_T_tau(tau[i]))*theta_d[i,:])

    Dc2_d[i,:] = c2eff_d[i,:]-ca2_tau(tau[i])

    # sigma_d[i,:] = ((1+1/3)*background_rho_nu_tau(tau[i])*one_time['shear_ur'])\
    #          /(background_rho_cdm_tau(tau[i])+(1+1/3)*background_rho_nu_tau(tau[i]))
print("Saving the files\n")

# for i in range(tau_num):
#     one_time = M_pre.get_transfer(background_z_at_tau(tau[i])) # transfer functions at each time tau

#     delta_d[i,:] = (d_ur[i]*background_rho_nu_tau(tau[i]) + background_rho_cdm_tau(tau[i])*d_cdm[i])\
#     /(background_rho_cdm_tau(tau[i])+background_rho_nu_tau(tau[i]))

# #     delta_T_std[i]= (background_rho_cdm_at_a_std(a_std[i])*delta_cdm[i]+background_rho_nu_at_a_std(a_std[i])*delta_nu[i])\
# #     /(background_rho_cdm_at_a_std(a_std[i])+background_rho_nu_at_a_std(a_std[i]))
# #     theta_T_std[i] = ((1+1/3)*background_rho_nu_at_a_std(a_std[i])*theta_nu[i])\
# #                     /((background_rho_cdm_at_a_std(a_std[i])+background_rho_nu_at_a_std(a_std[i])*(1.+1./3.)))

#     theta_d[i,:] = (1+1/3)*one_time['t_ur'][:]*background_rho_nu_tau(tau[i])\
#     /(background_rho_cdm_tau(tau[i])+(1+1/3)*background_rho_nu_tau(tau[i]))

#     # cs2_d[i,:] = (1./3.*d_ur[i]*background_rho_nu_tau(tau[i]))\
#     #         /(background_rho_cdm_tau(tau[i])*d_cdm[i]+background_rho_nu_tau(tau[i])*d_ur[i])

#     c2eff_d[i,:] = ((k*h)**2*cs2_d[i,:]*delta_d[i,:]+3*background_aH_at_tau(tau[i])*ca2_tau(tau[i])*(1+w_T_tau(tau[i]))*theta_d[i,:])\
#         /((k*h)**2*delta_d[i,:]+3*background_aH_at_tau(tau[i])*(1+w_T_tau(tau[i]))*theta_d[i,:])

#     Dc2_d[i,:] = c2eff_d[i,:]-ca2_tau(tau[i])

#     # sigma_d[i,:] = ((1+1/3)*background_rho_nu_tau(tau[i])*one_time['shear_ur'])\
#     #          /(background_rho_cdm_tau(tau[i])+(1+1/3)*background_rho_nu_tau(tau[i]))
# print("Saving the files\n")

np.savetxt(os.path.join(gdmdir,'k_values.dat'),k)
np.savetxt(os.path.join(gdmdir,'cs2_T.dat'), c2eff_d.ravel())
# np.savetxt(os.path.join(gdmdir,'sigma_T.dat'),sigma_d.ravel())
np.savetxt(os.path.join(gdmdir,'w_T.dat'), np.column_stack([a_tau,w_T_at_a_list,w_T_at_a_prime_list]))
