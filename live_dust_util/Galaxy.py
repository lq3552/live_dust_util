from .SnapshotContainer import SnapshotContainer
from .GrainSizeDistribution import GrainSizeDistribution
from .ExtinctionLaw import ExtinctionLaw
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
				   "GrainSizeDistributions",
				   "ExtinctionQuantities",
				   "VirialQuantities"]

	def __init__(self, snap):
		self.dataset = dict.fromkeys(Galaxy._field_list, None)
		self.dataset["MassesByType"] = dict.fromkeys(["PartType0", "PartType3", "PartType4"], 0.)
		for key in self.dataset["MassesByType"].keys():
			self.dataset["MassesByType"][key] = np.sum(snap.dataset[key + "/Masses"])
		self.dataset["Metallicity"] = np.sum(snap.dataset["PartType0/GFM_Metallicity"] * snap.dataset["PartType0/Masses"]) / self.dataset["MassesByType"]["PartType0"]
		self.dataset["GrainSizeDistributions"] = GrainSizeDistribution(snap)
		### self.dataset["ExtinctionQuantities"] ### I leave it to future work since I really don't like current implementation of ExtinctionLaw
			
	def show_field_list(self):
		print(list(self.dataset.keys()))
