import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.GrainSizeDistribution as GSD
import live_dust_util.ExtinctionLaw as Ext
import live_dust_util.ExtinctionLawParam as ExtPar

if __name__ == "__main__":
	"""
	test ExtinctionLaw and ExtinctionLawParam class, plot extinction curves
	"""
	gsd_tab = []

	a = 10**np.linspace(-3, 0 , 16)
	p_c = np.array([300, 300, 300])
	r_s = 0.0
	r_e = 30.0

	for i in range(30,36):
		snap = Snap(i,"MW_hr")
		gsd_tab.append(GSD(snap, a = a, p_c = p_c, r_s = r_s, r_e = r_e))

	wave = 10**np.linspace(-1,0,200)
	for gsd in gsd_tab:
		ext = Ext(gsd, wave, 
				  op_a = '../../grain_size.txt',
				  op_gra ='../../Gra_Optical/Gra_LD93_',
				  op_sil = '../../Sil_Optical/Sil_LD93_')
		plt.semilogx(wave, ext.extinction)
	
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=2.0),'k--',label='Cardelli')
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=3.1),'k--')
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=5.0),'k--')
	plt.semilogx(wave,ExtPar.smc(wave, tau_v=1), 'r--', label='Pei')
	plt.xlabel(r'$\lambda\ ({\rm \mu m})$')
	plt.ylabel(r'$A(\lambda)/A(V)$')
	plt.savefig('ext_mw_hr.png')
	plt.close()
