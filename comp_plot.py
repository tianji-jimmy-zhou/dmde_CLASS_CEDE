# import necessary modules
# uncomment to get plots displayed in notebook
%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from classy import Class
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import math

# esthetic definitions for the plots
font = {'size'   : 16, 'family':'STIXGeneral'}
axislabelfontsize='large'
matplotlib.rc('font', **font)
matplotlib.mathtext.rcParams['legend.fontsize']='medium'
plt.rcParams["figure.figsize"] = [8.0,6.0]

common_settings = {# which output? transfer functions only
                   'output':'tCl,pCl,lCl',
                   'l_max_scalars':3000,
                   'lensing':'yes',
                   # LambdaCDM parameters
                   'h':0.67556,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925,
                   'frac_dmde': 1e-10,
                   'g_scf': 0.0,
                   'm_scf': 0.0,
                   'Omega_k': 0,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.246,
                   # other output and precision parameters
                   'gauge':'synchronous'}

###############
#
# call CLASS
#
###############
M_lcdm = Class()
M_lcdm.set(common_settings)


M_lcdm.set({
'recfast_Nz0':100000,
'tol_thermo_integration':1.e-5,
'recfast_x_He0_trigger_delta':0.01,
'recfast_x_H0_trigger_delta':0.01,
'evolver':0,
'k_min_tau0':0.002,
'k_max_tau0_over_l_max':3.,
'k_step_sub':0.015,
'k_step_super':0.0001,
'k_step_super_reduction':0.1,
'start_small_k_at_tau_c_over_tau_h':0.0004,
'start_large_k_at_tau_h_over_tau_k':0.05,
'tight_coupling_trigger_tau_c_over_tau_h':0.005,
'tight_coupling_trigger_tau_c_over_tau_k':0.008,
'start_sources_at_tau_c_over_tau_h':0.006,
'l_max_g':50,
'l_max_pol_g':25,
'l_max_ur':50,
'tol_perturb_integration':1.e-6,
'perturb_sampling_stepsize':0.01,
'radiation_streaming_approximation':2,
'radiation_streaming_trigger_tau_over_tau_k':240.,
'radiation_streaming_trigger_tau_c_over_tau':100.,
'ur_fluid_approximation':2,
'ur_fluid_trigger_tau_over_tau_k':50.,
'l_logstep':1.026,
'l_linstep':25,
'hyper_sampling_flat':12.,
'hyper_nu_sampling_step':10.,
'hyper_phi_min_abs':1.e-10,
'hyper_x_tol':1.e-4,
'hyper_flat_approximation_nu':1.e6,
'q_linstep':0.20,
'q_logstep_spline':20.,
'q_logstep_trapzd':0.5,
'q_numstep_transition':250,
'transfer_neglect_delta_k_T_t2':100.,
'transfer_neglect_delta_k_T_e':100.,
'transfer_neglect_delta_k_T_b':100.,
'neglect_CMB_sources_below_visibility':1.e-30,
'transfer_neglect_late_source':3000.
})

M_lcdm.compute()

background = M_lcdm.get_background()

lcdm_z = background['z'] # read redshift

lcdm_H = background['H [1/Mpc]'] # read H
lcdm_rho_c = background['(.)rho_cdm'] # read rho_cdm
lcdm_DA = background['ang.diam.dist.'] # read DA
lcdm_rs = background['comov.snd.hrz.'] # read rs

lcdm_H_at_z = interp1d(lcdm_z,lcdm_H)
lcdm_rho_c_at_z = interp1d(lcdm_z,lcdm_rho_c)
lcdm_DA_at_z = interp1d(lcdm_z,lcdm_DA)
lcdm_rs_at_z = interp1d(lcdm_z,lcdm_rs)

lcdm_cl_lensed = M_lcdm.lensed_cl(3000)

M_lcdm.struct_cleanup()
M_lcdm.empty()

common_settings = {# which output? transfer functions only
                   'output':'tCl,pCl,lCl',
                   'l_max_scalars':3000,
                   'lensing':'yes',
                   # LambdaCDM parameters
                   'h':0.67556,
                   'omega_b':0.022032,
                   'omega_cdm':0.12038,
                   'A_s':2.215e-9,
                   'n_s':0.9619,
                   'tau_reio':0.0925,
                   'frac_dmde': 0.1,
                   'g_scf': 0.2,
                   'm_scf': 0.01,
                   'Omega_k': 0,
                   # Take fixed value for primordial Helium (instead of automatic BBN adjustment)
                   'YHe':0.246,
                   'gauge':'synchronous'}

###############
#
# call CLASS
#
###############
M_dmde = Class()
M_dmde.set(common_settings)

M_dmde.set({
'recfast_Nz0':100000,
'tol_thermo_integration':1.e-5,
'recfast_x_He0_trigger_delta':0.01,
'recfast_x_H0_trigger_delta':0.01,
'evolver':0,
'k_min_tau0':0.002,
'k_max_tau0_over_l_max':3.,
'k_step_sub':0.015,
'k_step_super':0.0001,
'k_step_super_reduction':0.1,
'start_small_k_at_tau_c_over_tau_h':0.0004,
'start_large_k_at_tau_h_over_tau_k':0.05,
'tight_coupling_trigger_tau_c_over_tau_h':0.005,
'tight_coupling_trigger_tau_c_over_tau_k':0.008,
'start_sources_at_tau_c_over_tau_h':0.006,
'l_max_g':50,
'l_max_pol_g':25,
'l_max_ur':50,
'tol_perturb_integration':1.e-6,
'perturb_sampling_stepsize':0.01,
'radiation_streaming_approximation':2,
'radiation_streaming_trigger_tau_over_tau_k':240.,
'radiation_streaming_trigger_tau_c_over_tau':100.,
'ur_fluid_approximation':2,
'ur_fluid_trigger_tau_over_tau_k':50.,
'l_logstep':1.026,
'l_linstep':25,
'hyper_sampling_flat':12.,
'hyper_nu_sampling_step':10.,
'hyper_phi_min_abs':1.e-10,
'hyper_x_tol':1.e-4,
'hyper_flat_approximation_nu':1.e6,
'q_linstep':0.20,
'q_logstep_spline':20.,
'q_logstep_trapzd':0.5,
'q_numstep_transition':250,
'transfer_neglect_delta_k_T_t2':100.,
'transfer_neglect_delta_k_T_e':100.,
'transfer_neglect_delta_k_T_b':100.,
'neglect_CMB_sources_below_visibility':1.e-30,
'transfer_neglect_late_source':3000.
})

M_dmde.compute()

background = M_dmde.get_background()

dmde_z = background['z'] # read redshift

dmde_H = background['H [1/Mpc]'] # read H
dmde_rho_c = background['(.)rho_cdm'] # read rho_cdm
dmde_DA = background['ang.diam.dist.'] # read DA
dmde_rs = background['comov.snd.hrz.'] # read rs

dmde_scf = background['phi_scf'] # read phi_scf

dmde_H_at_z = interp1d(dmde_z,dmde_H)
dmde_rho_c_at_z = interp1d(dmde_z,dmde_rho_c)
dmde_DA_at_z = interp1d(dmde_z,dmde_DA)
dmde_rs_at_z = interp1d(dmde_z,dmde_rs)
dmde_scf_at_z = interp1d(dmde_z,dmde_scf)

dmde_cl_lensed = M_dmde.lensed_cl(3000)

M_dmde.struct_cleanup()
M_dmde.empty()

plt.xlim([2,3000])

plt.xlabel(r"$\ell$")
plt.ylabel(r"$\Delta C_l^{XX}/C_l^{XX,{\Lambda{\rm CDM}}}$")
plt.title(r"$r=0.1$")
plt.grid()

ell = lcdm_cl_lensed['ell']
factor = 1.e10*ell*(ell+1.)/2./math.pi

plt.semilogx(ell,(lcdm_cl_lensed['tt']-dmde_cl_lensed['tt'])/lcdm_cl_lensed['tt'],'r-',label=r'$\mathrm{TT}$')
plt.semilogx(ell,(lcdm_cl_lensed['ee']-dmde_cl_lensed['ee'])/lcdm_cl_lensed['ee'],'b:',label=r'$\mathrm{EE}$')

plt.legend(loc='right',bbox_to_anchor=(1.4, 0.5))

plt.savefig('cl_comp.pdf',bbox_inches='tight')
