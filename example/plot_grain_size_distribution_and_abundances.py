import numpy as np
import matplotlib.pyplot as plt
import live_dust_util.SnapshotContainer as Snap
import live_dust_util.GrainSizeDistribution as GSD

if __name__ == "__main__":
	"""
	test GrainSizeDistribution class
	"""

	"""
	plot grain size distributions
	"""
	gsd_tab = []
	for i in range(50,56):
		snap = Snap(i,"SMC_hr")
		gsd_tab.append(GSD(snap))

	i = 0
	for gsd in gsd_tab:
	
		fig, ax = plt.subplots(2,1, sharex = True)
	
		for key in GSD.species_keys:
			ax[0].loglog(gsd.a, gsd.DMSF[key], '-', label=key)
		ax[0].legend()
		ax[0].set_ylabel(r'$\rho_{\rm gr} a^3 \times N $')
	
		ax[1].loglog(gsd.a,gsd.DMSF['PAH'] / (gsd.DMSF['Aliphatic C'] + gsd.DMSF['PAH']) ,'-',)
		ax[1].set_xlabel(r'$a\ (\mu m)$')
		ax[1].set_ylabel(r'PAH / Carbonaceous')
		ax[1].set_ylim(0.1, 1.1)
	
		plt.savefig(f'dmsf_smc_{i: 03d}.png',dpi=300)
		plt.close()
		i += 1
	
	
	"""
	plot evolution of grain species abundances
	"""
	gsd_tab = []
	for i in range(6,51):
		snap = Snap(i,"SMC_hr")
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
