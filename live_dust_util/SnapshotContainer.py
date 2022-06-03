import h5py

class SnapshotContainer(object):
	"""
	Container of some important fields extracted from Arepo_dust simulation snapshots
	To show list of fields:
		use show_field_list()

	Parameter:
		snap_no: <int> No. of the input snapshot
		snap_dir: <str> directory name of the input snapshot; default = '.'
		snap_pref: <str> prefix of the input snapshot; default = 'snapshot'
	"""
	__field_list = ["PartType0/Density",
	                "PartType0/Masses",
					"PartType0/Coordinates",
					"PartType0/GFM_Metallicity",
					"PartType0/GFM_Metals",
					"PartType0/MolecularHFrac",
					"PartType0/StarFormationRate",
		            "PartType3/Dust_NumGrains",
		            "PartType3/Dust_DustDensity",
					"PartType3/Masses",
					"PartType3/Coordinates",
					"PartType3/Dust_MetalFractions",
					"PartType4/Masses",
					"PartType4/SNIaNumber",
					"PartType4/SNIINumber"] # I seriously don't want users to change it

	def __init__(self, snap_no, snap_dir = '.', snap_pref = 'snapshot'):
		self.dataset = dict.fromkeys(SnapshotContainer.__field_list, None)
		snap = h5py.File('%s/%s_%03d.hdf5' % (snap_dir, snap_pref, snap_no), 'r')
		for key in self.dataset.keys():
			self.dataset[key] = snap[key][()]
			
	def show_field_list(self):
		print(list(self.dataset.keys()))
