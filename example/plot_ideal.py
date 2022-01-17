import sys
import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.GrainSizeDistribution as GSD
import live_dust_util.ExtinctionLaw as Ext
import live_dust_util.ExtinctionLawParam as ExtPar

'''
gsd_tab = []
for i in range(6,50):
	snap = Snap(i,'../../smc_mr_fuv')
	gsd_tab.append(GSD(snap))
i = 6
for gsd in gsd_tab:
	print(i)

	fig, ax = plt.subplots(2,1, sharex = True)

	for key in GSD.species_keys:
		ax[0].loglog(gsd.a,gsd.DMSF[key],'-',label=key)
	ax[0].legend()
	ax[0].set_ylabel(r'$\rho_{\rm gr} a^3 \times N $')

	ax[1].semilogx(gsd.a,gsd.DMSF['PAH'] / (gsd.DMSF['Aliphatic C'] + gsd.DMSF['PAH']) ,'-',)
	ax[1].set_xlabel(r'$a\ (\mu m)$')
	ax[1].set_ylabel(r'PAH / Carbonaceous')
	ax[1].set_ylim(0.5, 1.1)

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


gsd_tab = []
for i in range(6,49):
	snap = Snap(i,'../../smc_mr_fuv')
	gsd_tab.append(GSD(snap))
C_frac, PAH_frac, Si_frac = [], [], []
tsim = []
i = 0
for gsd in gsd_tab:
	abu = gsd.compute_abundances()
	C_frac.append(abu['Aliphatic C'] + abu['PAH'])
	PAH_frac.append(abu['PAH'])
	Si_frac.append(abu['Silicate'])
	tsim.append(i * 0.05)
	i += 1
plt.plot(tsim,C_frac,label='Carbonaceous')
plt.plot(tsim,PAH_frac,label='PAH')
plt.plot(tsim,Si_frac,label='Silicate')

plt.xlabel(r'$t$ (Gyr)')
plt.ylabel(r'Abundances')
plt.legend()
plt.savefig('C_abu.png')
plt.close()

