import sys
import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.GrainSizeDistribution as GSD
import live_dust_util.ExtinctionLaw as Ext
import live_dust_util.ExtinctionLawParam as ExtPar

'''
gsd_tab = []
for i in range(6,200):
	snap = Snap(i,'../../smc_mr_fuv')
	gsd_tab.append(GSD(snap))
i = 6
for gsd in gsd_tab:
	print(i)
	for key in GSD.species_keys:
		plt.loglog(gsd.a,gsd.DMSF[key],'--',label=key)
		
	plt.legend()
	plt.xlabel(r'$a\ (\mu m)$')
	plt.ylabel(r'$\rho_{\rm gr} a^3 \times N $')
	plt.savefig('dmsf_smc.%03d.png' % i,dpi=300)
	plt.close()
	i += 1
'''

gsd_tab = []
for i in range(30,31):
	snap = Snap(i,'../../smc_mr_fuv')
	gsd_tab.append(GSD(snap))
#for i in range(37,38):
#	snap = Snap(i,'../../mw_lr_fuv')
#	gsd_tab.append(GSD(snap))
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

'''
snap_tab = []
gsd_tab = []
for i in range(6,200):
	snap = Snap(i,'../../smc_mr_fuv')
	snap_tab.append(snap)
	gsd_tab.append(GSD(snap))
Cfrac = []
tsim = []
i = 0
for gsd in gsd_tab:
	snap = snap_tab[i]
	abu = gsd.compute_abundances()
	Cfrac.append(abu['Aliphatic C'] + abu['PAH'])
	tsim.append(i * 0.05)
	i += 1
plt.plot(tsim,Cfrac,label='SMC_PAH')

plt.xlabel(r'$t$ (Gyr)')
plt.ylabel(r'Carbonaceous / Total')
plt.legend()
plt.savefig('C_abu.png')
plt.close()
'''
