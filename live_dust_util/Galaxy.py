from .SnapshotContainer import SnapshotContainer
from .GrainSizeDistribution import GrainSizeDistribution
from .ExtinctionLaw import ExtinctionLaw
from .utils import IndexUtil as iutil
from .utils import MoleUtil as mutil
import numpy as np

class Galaxy(object):
	"""
	Container of galaxy informations
	To show list of fields:
		use show_field_list() method

	Parameter:
		snap:  <SnapshotContainer>
		[p_c]: <ndarray[3], dtype = float64> center to compute radii in code units
		[r_s]: <float32> lower bound of radius interval in code units
		[r_e]: <float32> upper bound of radius interval in code units
		[lz]:  <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
	"""
	_field_list = [
				   "MassesByType",
				   "Metallicity",
				   "SFR",
				   "RatioCOtoC",
				   "GrainSizeDistributions",
				   "ExtinctionQuantities",
				   "VirialQuantities",
				   "SNNumber",
				   "TauShatter",
				   "TauCoagulate"
				  ]

	def __init__(self, snap, p_c = [], r_s = None, r_e = None, lz = np.array([0.,0.,1.])):
		self.filt = snap.compute_filter(p_c, r_s, r_e, lz)
		self.dataset = dict.fromkeys(Galaxy._field_list, None)
		self.dataset["MassesByType"] = dict.fromkeys(["PartType0", "PartType3", "PartType4"], 0.)
		for key in self.dataset["MassesByType"].keys():
			self.dataset["MassesByType"][key] = np.sum(snap.dataset[key + "/Masses"][self.filt[key]])
		self.dataset["Metallicity"] = np.sum(snap.dataset["PartType0/GFM_Metallicity"][self.filt["PartType0"]]
				                           * snap.dataset["PartType0/Masses"][self.filt["PartType0"]]) \
										   / self.dataset["MassesByType"]["PartType0"]
		self.dataset["SNNumber"] = np.sum(snap.dataset["PartType4/SNIaNumber"][self.filt["PartType4"]]
								        + snap.dataset["PartType4/SNIINumber"][self.filt["PartType4"]])
		self.dataset["SFR"] = np.sum(snap.dataset["PartType0/StarFormationRate"][self.filt["PartType0"]]) 
		tauShatter = snap.dataset["PartType3/Dust_TauShatter"][self.filt["PartType3"]]
		mShatter = snap.dataset["PartType3/Masses"][self.filt["PartType3"]]
		filt2 = np.where(tauShatter < 1e2)
		tauShatter = tauShatter[filt2]
		mShatter = mShatter[filt2]
		tauCoagulate = snap.dataset["PartType3/Dust_TauCoagulate"][self.filt["PartType3"]]
		mCoagulate = snap.dataset["PartType3/Masses"][self.filt["PartType3"]]
		filt2 = np.where(tauCoagulate < 1e2)
		tauCoagulate = tauCoagulate[filt2]
		mCoagulate = mCoagulate[filt2]
		self.dataset["TauShatter"] = np.sum(tauShatter * mShatter) / np.sum(mShatter)
		self.dataset["TauCoagulate"] =  np.sum(tauCoagulate * mCoagulate) / np.sum(mCoagulate)
		### self.dataset["GrainSizeDistributions"] = GrainSizeDistribution(snap) ###
		### self.dataset["ExtinctionQuantities"] ### I leave it to future work since I really don't like current implementation of ExtinctionLaw
		self._compute_abundances(snap)
			
	def show_field_list(self):
		print(list(self.dataset.keys()))

	def _compute_abundances(self, snap):
		NC = snap.dataset["PartType0/GFM_Metals"][self.filt["PartType0"], iutil.elem_i_C] / mutil.A_C
		NO = snap.dataset["PartType0/GFM_Metals"][self.filt["PartType0"], iutil.elem_i_O] / mutil.A_O
		self.dataset["RatioCOtoC"] = np.sum(np.minimum(NC, NO) * snap.dataset["PartType0/MolecularHFrac"][self.filt["PartType0"]] * 0.5 * mutil.A_C)\
			                         / np.sum(snap.dataset["PartType0/GFM_Metals"][self.filt["PartType0"], iutil.elem_i_C])
