import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.GrainSizeDistribution as GSD
import live_dust_util.ExtinctionLaw as Ext
import live_dust_util.ExtinctionLawParam as ExtPar

if __name__ == "__main__":
	"""
	test ExtinctionLaw class, plot extinction curves
	"""
	gsd_tab = []
	for i in range(30,31):
		snap = Snap(i,'../../smc_mr_h50pc')
		gsd_tab.append(GSD(snap))
	for i in range(20,21):
		snap = Snap(i,'../../mw_mr_h100pc')
		gsd_tab.append(GSD(snap))
	wave = 10**np.linspace(-1,0,200)
	label=['SMC','MW']
	i = 0
	for gsd in gsd_tab:
		ext = Ext(gsd, wave, 
				  op_a = '../../grain_size.txt',
				  op_gra ='../../Gra_Optical/Gra_LD93_',
				  op_sil = '../../Sil_Optical/Sil_LD93_')
		plt.semilogx(wave, ext.extinction, label=label[i])
		i += 1
	
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=2.0),'k--',label='Cardelli')
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=3.1),'k--')
	plt.semilogx(wave,ExtPar.cardelli(wave, tau_v=1, R_v=5.0),'k--')
	plt.semilogx(wave,ExtPar.smc(wave, tau_v=1), 'r--', label='Pei')
	plt.xlabel(r'$\lambda\ ({\rm \mu m})$')
	plt.ylabel(r'$A(\lambda)/A(V)$')
	plt.legend()
	plt.savefig('ext.png')
	plt.close()
