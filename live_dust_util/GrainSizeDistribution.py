import copy
import numpy as np
from .SnapshotContainer import SnapshotContainer
from .utils import IndexUtil

class GrainSizeDistribution(object):
	"""
	Contains the grain number count or mass as a function of grain sizes in one galaxy.

	Parameters:
		snap: <SnapshotContainer>
		a  : <ndarray[], dtype = float64> centers of grain size bins
		p_c: <ndarray[3], dtype = float64> center to compute radii
		r_s: <float32> lower bound of radius interval
		r_e: <float32> upper bound of radius interval
		lz: <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
	"""
	species_rho  = {"Aliphatic C": 2.2, "PAH": 2.2, "Silicate": 3.3}
	species_keys = species_rho.keys()


	def __init__(self, snap, a = 10**np.linspace(-3,0,16), p_c = [], r_s = None, r_e = None, lz = np.array([0.,0.,1.])):
		self.set_grain_size_distribution(snap, a, p_c, r_s, r_e, lz)


	def set_grain_size_distribution(self,snap, a, p_c, r_s, r_e, lz):
		"""
		Computes properties related to grain size distributions from a snapshot

		Parameters:
			snap: <SnapshotContainer>
			a  : <ndarray[], dtype = float64> centers of grain size bins
			p_c: <ndarray[3], dtype = float64> center to compute radii
			r_s: <float32> lower bound of radius interval
			r_e: <float32> upper bound of radius interval
			lz: <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
		"""

		self.a = a # I need to refactor this because ideally 
		if len(self.a) > 1:
			self.dloga = np.log10(self.a[1] / self.a[0])
#		the method should be able to read grain_size_distributions directly from the 
#		snapshot (which should be saved as a subdataset of PartType3 or an attribute 
#		of the header but as a practice I simply ignore this for now. So does number
#		of species.
		self.DNSF = dict.fromkeys(GrainSizeDistribution.species_keys, None)
		self.DMSF = copy.deepcopy(self.DNSF)

		nbins = len(self.a)

		if (p_c == []) or (r_s == None) or (r_e == None):
			filt = np.where(snap.dataset['PartType3/Masses'] > 0) # essentially no filter applied
		else:
			x = snap.dataset["PartType3/Coordinates"][:,0] - p_c[0]
			y = snap.dataset["PartType3/Coordinates"][:,1] - p_c[1]
			z = snap.dataset["PartType3/Coordinates"][:,2] - p_c[2]
			r_s2 = r_s**2
			r_e2 = r_e**2
			r2 = x**2 + y**2 + z**2
			filt = np.where((r2 >= r_s2) & (r2 < r_e2))
		f_PAH = snap.dataset["PartType3/Dust_NumGrains"][filt][:,-nbins:]
		mass_silicate_test = snap.dataset["PartType3/Dust_MetalFractions"][filt][:, IndexUtil.elem_i_Si] / 0.166
		mass_silicate_mix  = 1.0 - snap.dataset["PartType3/Dust_MetalFractions"][filt][:, IndexUtil.elem_i_C];
		self.DNSF["Aliphatic C"] = np.sum((1.0 - f_PAH)
			                       * snap.dataset["PartType3/Dust_NumGrains"][filt][:,nbins: 2 * nbins], axis=0)
		self.DNSF["PAH"] = np.sum(f_PAH * snap.dataset["PartType3/Dust_NumGrains"][filt][:,nbins: 2 * nbins], axis=0)
		self.DNSF["Silicate"] = np.matmul(mass_silicate_test / mass_silicate_mix,\
					            snap.dataset["PartType3/Dust_NumGrains"][filt][:, : nbins])
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
		return arr * 4 * np.pi / 3 * self.a**3 * GrainSizeDistribution.species_rho[key] # cgs
