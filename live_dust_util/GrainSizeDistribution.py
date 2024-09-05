import copy
import numpy as np
from .SnapshotContainer import SnapshotContainer
from .utils import IndexUtil

class GrainSizeDistribution(object):
	"""
	Contains the grain number count or mass as a function of grain sizes in one galaxy.

	Parameters:
		snap:  <SnapshotContainer>
		[a]:   <ndarray[], dtype = float64> centers of grain size bins in micron
		[p_c]: <ndarray[3], dtype = float64> center to compute radii in code units
		[r_s]: <float32> lower bound of radius interval in code units
		[r_e]: <float32> upper bound of radius interval in code units
		[lz]:  <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
	"""
	# class variables
	species_rho      = {"Aliphatic C": 2.2, "PAH": 2.2, "Silicate": 3.3}
	species_rho_ext  = {"Aliphatic C": 2.2, "PAH": 2.2, "Silicate": 3.3, "Carbonaceous": 3.3} #better solution to this?
	species_keys     = species_rho.keys()
	species_keys_ext = species_rho_ext.keys()
	DL_GRAIN_BINS = 16
	SEC_PER_GIGAYEAR = 3.15576e16
	PROTONMASS = 1.67262178e-24
	CM_TO_UM = 1e4


	def __init__(self, snap, a = 10**np.linspace(-3,0,16), p_c = [], r_s = None, r_e = None, lz = np.array([0.,0.,1.])):
		self.set_grain_size_distribution(snap, a, p_c, r_s, r_e, lz)


	def set_grain_size_distribution(self,snap, a, p_c, r_s, r_e, lz):
		"""
		Computes properties related to grain size distributions from a snapshot

		Parameters:
			snap: <SnapshotContainer>
			a  : <ndarray[], dtype = float64> centers of grain size bins in micron
			p_c: <ndarray[3], dtype = float64> center to compute radii in code units
			r_s: <float32> lower bound of radius interval in code units
			r_e: <float32> upper bound of radius interval in code units
			lz:  <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
		"""
		self.snap = snap
		self.a = a # I need to refactor this because ideally 
		if len(self.a) > 1:
			self.dloga = np.log10(self.a[1] / self.a[0])
#		the method should be able to read grain_size_distributions directly from the 
#		snapshot (which should be saved as a subdataset of PartType3 or an attribute 
#		of the header but as a practice I simply ignore this for now. So does number
#		of species.
		self.DNSF = dict.fromkeys(GrainSizeDistribution.species_keys_ext, None)
		self.DMSF = copy.deepcopy(self.DNSF)

		nbins = len(self.a)

		filt = snap.compute_filter(p_c, r_s, r_e, lz)['PartType3']
		if (snap.dataset["PartType3/Dust_NumGrains"].shape[1]) >= 3 * nbins:
			f_PAH = snap.dataset["PartType3/Dust_NumGrains"][filt][:,-nbins:]
			self.DNSF["Aliphatic C"] = np.sum((1.0 - f_PAH)
									  * snap.dataset["PartType3/Dust_NumGrains"][filt][:,nbins: 2 * nbins], axis=0)
			self.DNSF["PAH"] = np.sum(f_PAH * snap.dataset["PartType3/Dust_NumGrains"][filt][:,nbins: 2 * nbins], axis=0)
			self.DNSF["Carbonaceous"] = np.sum(snap.dataset["PartType3/Dust_NumGrains"][filt][:,nbins: 2 * nbins], axis=0)
			self.DNSF["Silicate"] = np.sum(snap.dataset["PartType3/Dust_NumGrains"][filt][:, : nbins], axis=0)
		else:
			f_C = snap.dataset["PartType3/Dust_MetalFractions"][filt][: ,IndexUtil.elem_i_C]
			self.DNSF["Aliphatic C"] = np.sum(np.dot(f_C, snap.dataset["PartType3/Dust_NumGrains"][filt][:, :]), axis=0)
			self.DNSF["Silicate"] = np.sum(np.dot((1.0 - f_C), snap.dataset["PartType3/Dust_NumGrains"][filt][:, :]), axis=0)
			self.DNSF["PAH"] = np.sum(0.0 * snap.dataset["PartType3/Dust_NumGrains"][filt][:, :], axis=0)
			self.DNSF["Carbonaceous"] = self.DNSF["Aliphatic C"] # no need to deep copy since they are actually the same
		for key in GrainSizeDistribution.species_keys_ext:
			self.DMSF[key] = self._from_n_to_m(self.DNSF[key],key)
		'''
		# compute rates of physical processes
		self.rates = dict.fromkeys(SnapshotContainer.__rates_list, np.zeros((snap.dataset["PartType3/Masses"].shape[0], 2))) ### TODO: more flexible
		self._compute_rates(0)
		self._compute_rates(1)
		'''

	def get_grain_size_distribution(self, spe, qtype = "mass"):
		"""
		return the grain number count or mass as a function of grain sizes.

		Parameters:
		spe: <str> species "Aliphatic C", "PAH" or "Silicate"
		qtype: <str> "mass" or "num", default = "mass"
		"""
		if spe in GrainSizeDistribution.species_keys_ext:
			if qtype in ["Mass", "mass"]:
				return self.DMSF[spe]
			elif qtype in ["Num", "num"]:
				return self.DNSF[spe]
			else:
				print("Field type not found! Please make sure:")
				print("field type keyword in", ["mass", "num"])
		else:
			print("Species not found! Please make sure:")
			print("species keyword in",list(GrainSizeDistribution.species_keys_ext))
			return np.array([])


	def compute_small_to_large_ratio(self):
		"""
		compute small-to-large-grain mass ratio for different grain species

		return: <dict> small-to-large mass ratio
		"""
		i = 0
		stl = np.zeros(len(GrainSizeDistribution.species_keys_ext))
		for key in GrainSizeDistribution.species_keys_ext:
			filt_small = np.where(self.a <= 6e-2) # Aoyama+2020
			filt_large = np.where(self.a > 6e-2)
			m_small = np.sum(self.DMSF[key][filt_small])
			m_large = np.sum(self.DMSF[key][filt_large])
			stl[i] = m_small / m_large
			i += 1
		return dict(zip(GrainSizeDistribution.species_keys_ext, stl))

	def compute_abundances(self):
		"""
		compute abundances of different grain species

		return: <dict> abundances
		"""
		m_spe = np.zeros(len(GrainSizeDistribution.species_keys_ext))
		i = 0
		for key in GrainSizeDistribution.species_keys_ext:
			m_spe[i] = np.sum(self.DMSF[key])
			i += 1
		return dict(zip(GrainSizeDistribution.species_keys_ext, m_spe / np.sum(m_spe[:-1])))

	def _from_n_to_m(self, arr, key):
		return arr * 4 * np.pi / 3 * self.a**3 * GrainSizeDistribution.species_rho_ext[key] # cgs

	'''
	def _compute_rates(self, s):
		self.__compute_rates_shatter_coagulate(s)
		#self.rates["Shattering"] = 0
		#self.rates["Coagulation"] = 0
		print(self.rates)


	def __compute_rates_shatter_coagulate(self, s)
		Mass = self.snap.dataset["PartType3/Masses"]
		DustDensity = self.snap.dataset["PartType3/Dust_DustDensity"]
		cf_atime = 1 / (1 + self.snap.header["Redshift"])
		cf_a3inv = 1 / cf_atime ** 3
		print(cf_a3inv)
		UnitLength_in_cm = self.snap.header["UnitLength_in_cm"]
		UnitMass_in_g = self.snap.header["UnitMass_in_g"]
		UnitDensity_in_cgs = UnitMass_in_g / UnitLength_in_cm ** 3
		UnitVelocity_in_cm_per_s = self.snap.header["UnitVelocity_in_cm_per_s"]
		HubbleParam = self.snap.header["HubbleParam"]

		V_d = Mass / (DustDensity * cf_a3inv)
		V_d_um3 = V_d * pow(UnitLength_in_cm / HubbleParam * self.CM_TO_UM, 3.0)

		v_rels = np.zeros((self.DL_GRAIN_BINS, self.DL_GRAIN_BINS))
		v_grain = np.zeros(self.DL_GRAIN_BINS)
		cs = self.snap.dataset["PartType3/Dust_GasSoundSpeed"] * np.sqrt(All.cf_atime) * UnitVelocity_in_cm_per_s
		T_local = self.snap.dataset["PartType3/Dust_GasTemperature"]
		rho_local = self.snap.dataset["PartType3/Dust_GasDensity"] * cf_a3inv
		rho_ref = self.PROTONMASS / (UnitDensity_in_cgs * HubbleParam * HubbleParam)
		
	'''	
