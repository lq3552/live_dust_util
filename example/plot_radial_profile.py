import sys
import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.RadialProfile as RadP

hubble = 0.7
Lcode = 3.085678e21
Mcode = 1.989e43 
Msun = 1.989e33
Vcode = 1e5
Dcode = (Mcode / hubble) / (Lcode / hubble)**3
Mp = 1.674e-24

a = 10**np.linspace(-3,0,16)
p_c = np.array([300, 300, 300])
r_s = 0.0
r_e = 30.0
n_bins = 14
wave = 10**np.linspace(-1,0,200)

fig, ax = plt.subplots(nrows = 4, ncols = 2, sharex = True, figsize=(5*2, 4*4))

snap = Snap(40, "./MW_hr")
rad_prof = RadP(snap, wave,
					op_a = './grain_size.txt',
					op_gra ='./Gra_Optical/Gra_LD93_',
					op_sil = './Sil_Optical/Sil_LD93_',
					a = a,
					p_c = p_c,
					r_s = r_s,
					r_e = r_e,
					n_bins = n_bins
				)
ax[0,0].semilogy(rad_prof.rad, rad_prof.data_radial["GasDensity"] * Dcode / Mp, 'o', label = "Photoionization Feedback") #label = r"$Y_{\rm C}$ $\times 2$")
ax[1,0].semilogy(rad_prof.rad, rad_prof.data_radial["DustDensity"] * Dcode / Mp , 'o')
ax[2,0].semilogy(rad_prof.rad, rad_prof.data_radial["DustDensity"]/rad_prof.data_radial["GasDensity"], 'o')
ax[3,0].semilogy(rad_prof.rad, rad_prof.data_radial["SigmaGas"] / rad_prof.data_radial["SigmaSFR"], 'o')
ax[0,1].plot(rad_prof.rad, rad_prof.data_radial["BumpStrength"], 'o')
ax[1,1].plot(rad_prof.rad, rad_prof.data_radial["UVOPSlope"], 'o')
ax[2,1].plot(rad_prof.rad, rad_prof.data_radial["RV"], 'o')
ax[3,1].plot(rad_prof.rad, rad_prof.data_radial["CarbonFraction"], 'o')

ax[0,0].set_ylabel(r"$\rho_{\rm g}$ ($m_{\rm p}\ {\rm cm}^{-3}$)")
ax[0,0].legend()
ax[1,0].set_ylabel(r"$\rho_{\rm d}$ ($m_{\rm p}\ {\rm cm}^{-3}$)")
ax[2,0].set_ylabel("DTG")
ax[3,0].set_ylabel(r"$\tau_{\rm SFR}$ (${\rm yr}$)")
ax[3,0].set_xlabel(r"$r$ (kpc)")
ax[0,1].set_ylabel("Bump")
ax[1,1].set_ylabel("Slope")
ax[2,1].set_ylabel(r"$R_V$")
ax[3,1].set_ylabel(r"$f_{\rm C}$")
ax[3,1].set_xlabel(r"$r$ (kpc)")
plt.savefig("radial_mw_hr.png",dpi=300)
plt.close()
