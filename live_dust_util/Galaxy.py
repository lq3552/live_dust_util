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
		snap: <SnapshotContainer>
	"""
	_field_list = ["MassesByType",
				   "Metallicity",
				   "RatioCOtoC",
				   "GrainSizeDistributions",
				   "ExtinctionQuantities",
				   "VirialQuantities",
				   "SNNumber"]

	def __init__(self, snap):
		self.dataset = dict.fromkeys(Galaxy._field_list, None)
		self.dataset["MassesByType"] = dict.fromkeys(["PartType0", "PartType3", "PartType4"], 0.)
		for key in self.dataset["MassesByType"].keys():
			self.dataset["MassesByType"][key] = np.sum(snap.dataset[key + "/Masses"])
		self.dataset["Metallicity"] = np.sum(snap.dataset["PartType0/GFM_Metallicity"] 
				                      * snap.dataset["PartType0/Masses"]) / self.dataset["MassesByType"]["PartType0"]
		self.dataset["GrainSizeDistributions"] = GrainSizeDistribution(snap)
		self.dataset["SNNumber"] = np.sum(snap.dataset["PartType4/SNIaNumber"]+snap.dataset["PartType4/SNIINumber"])
		### self.dataset["ExtinctionQuantities"] ### I leave it to future work since I really don't like current implementation of ExtinctionLaw
		self._compute_abundances(snap)
			
	def show_field_list(self):
		print(list(self.dataset.keys()))

	def _compute_abundances(self, snap):
		NC = snap.dataset["PartType0/GFM_Metals"][:, iutil.elem_i_C] / mutil.A_C
		NO = snap.dataset["PartType0/GFM_Metals"][:, iutil.elem_i_O] / mutil.A_O
		self.dataset["RatioCOtoC"] = np.sum(np.minimum(NC, NO) * snap.dataset["PartType0/MolecularHFrac"] * 0.5 * mutil.A_C)\
			                         / np.sum(snap.dataset["PartType0/GFM_Metals"][:, iutil.elem_i_C])
