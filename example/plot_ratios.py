import live_dust_util.SnapshotContainer as Snap
from live_dust_util import Galaxy
from live_dust_util import GalaxyCatalog
import matplotlib.pyplot as plt
import numpy as np


catalog = GalaxyCatalog()
for i in range(4,39):
	snap = Snap(i,'../../mw_lr_fuv')
	gal  = Galaxy(snap)
	catalog.add(gal)


Z    = np.array([gal.dataset["Metallicity"] for gal in catalog]) / 0.0127
print(Z)
Md   = np.array([gal.dataset["MassesByType"]["PartType3"] for gal in catalog])
qPAH = np.array([gal.dataset["GrainSizeDistributions"].compute_abundances()["PAH"] for gal in catalog])
Mg   = np.array([gal.dataset["MassesByType"]["PartType0"] for gal in catalog])
Mpah = Md * qPAH

plt.loglog(Z, qPAH, 'o', label="qPAH")
plt.loglog(Z, Md/Mg, 'o', label="DTG")
plt.legend()
plt.xlabel(r"$Z/Z_{\odot}$")
plt.ylabel("ratios")
plt.show()
