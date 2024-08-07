import live_dust_util.SnapshotContainer as Snap
from live_dust_util import Galaxy
from live_dust_util import GalaxyCatalog
import matplotlib.pyplot as plt
import numpy as np


a = 10**np.linspace(-3, 0 , 16)
p_c = np.array([300, 300, 300])
r_s = 0.0
r_e = 5.0

catalog = GalaxyCatalog()
for i in range(6,59):
	snap = Snap(i,"SMC_hr")
	gal  = Galaxy(snap, a = a, p_c = p_c, r_s = r_s, r_e = r_e)
	catalog.add(gal)


Z    = np.array([gal.dataset["Metallicity"] for gal in catalog]) / 0.0127
print(Z)
Md   = np.array([gal.dataset["MassesByType"]["PartType3"] for gal in catalog])
qPAH = np.array([gal.dataset["GrainSizeDistributions"].compute_abundances()["PAH"] for gal in catalog])
Mg   = np.array([gal.dataset["MassesByType"]["PartType0"] for gal in catalog])
Mpah = Md * qPAH
fCO  = np.array([gal.dataset["RatioCOtoC"] for gal in catalog])

plt.loglog(Z, qPAH, 'o', label="qPAH")
plt.loglog(Z, Md/Mg, 'o', label="DTG")
plt.loglog(Z, fCO, 'o', label=r"f$_{\rm CO}$")
plt.legend()
plt.xlabel(r"$Z/Z_{\odot}$")
plt.ylabel("ratios")
plt.show()
