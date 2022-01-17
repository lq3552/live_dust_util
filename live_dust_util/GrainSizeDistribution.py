import copy
import numpy as np
from .SnapshotContainer import SnapshotContainer

class GrainSizeDistribution(object):
	"""
	Contains the grain number count or mass as a function of grain sizes in one galaxy.

	Parameters:
		snap: <SnapshotContainer>
	"""
	species_rho  = {"Aliphatic C": 2.2, "PAH": 2.2, "Silicate": 3.3}
	species_keys = species_rho.keys()


	def __init__(self, snap):
		self.set_grain_size_distribution(snap)


	def set_grain_size_distribution(self, snap):
		"""
		Computes properties related to grain size distributions from a snapshot

		Parameters:
			snap: <SnapshotContainer>
		"""

		self.a = 10** np.linspace(-3,0,16) # I need to refactor this because ideally 
#		the method should be able to read grain_size_distributions directly from the 
#		snapshot (which should be saved as a subdataset of PartType3 or an attribute 
#		of the header but as a practice I simply ignore this for now. So does number
#		of species.
		self.DNSF = dict.fromkeys(GrainSizeDistribution.species_keys, [])
		self.DMSF = copy.deepcopy(self.DNSF)

		nbins = len(self.a)

		f_PAH = snap.dataset['PartType3/Dust_NumGrains'][:,-nbins:]
		self.DNSF["Aliphatic C"] = np.sum((1.0 - f_PAH)
			                       * snap.dataset['PartType3/Dust_NumGrains'][:,nbins: 2 * nbins], axis=0)
		self.DNSF["PAH"] = np.sum(f_PAH * snap.dataset['PartType3/Dust_NumGrains'][:,nbins: 2 * nbins], axis=0)
		self.DNSF["Silicate"] = np.sum(snap.dataset['PartType3/Dust_NumGrains'][:, : nbins], axis=0)
		for key in GrainSizeDistribution.species_keys:
			self.DMSF[key] = self._from_n_to_m(self.DNSF[key],key)


	def get_grain_size_distribution(self, spe, qtype):
		"""
		return the grain number count or mass as a function of grain sizes.

		Parameters:
		spe: str, species: "Aliphatic C", "PAH" or "Silicate"
		qtype: "mass" or "num"
		"""
		if spe in GrainSizeDistribution.species_keys:
			if qtype in ["Mass", "mass"]:
				return self.DMSF[spe]
			elif qtype in ["Num", "num"]:
				return self.DNSF[spe]
			else:
				print("Field type not found! Please make sure:")
				print("field type keyword in", ["mass", "num"])
		else:
			print("Species not found! Please make sure:")
			print("species keyword in",list(GrainSizeDistribution.species_keys))
			return np.array([])


	def compute_abundances(self):
		"""
		compute abundances of different grain species

		return: dict, abundances
		"""
		n_spe = len(GrainSizeDistribution.species_keys)
		m_spe = np.zeros(n_spe)
		i = 0
		for key in GrainSizeDistribution.species_keys:
			m_spe[i] = np.sum(self.DMSF[key])
			i += 1

		return dict(zip(GrainSizeDistribution.species_keys,m_spe / np.sum(m_spe)))


	def _from_n_to_m(self, arr, key): # TODO: this looks ugly
		return arr * self.a**3 * GrainSizeDistribution.species_rho[key] # cgs
