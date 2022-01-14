import h5py

class SnapshotContainer(object):
	"""
	Container of some important fields extracted from Arepo_dust simulation snapshots
	List of fields:
		use show_field_list

	Parameter:
		snap_no: <int> No. of the input snapshot
		snap_dir: <str> directory name of the input snapshot; default = '.'
		snap_pref: <str> prefix of the input snapshot; default = 'snapshot'
	"""
	__field_list = ['PartType3/Dust_NumGrains',
					'PartType3/Masses',
					'PartType4/Masses'] # I seriously don't want you to change it

	def __init__(self, snap_no, snap_dir = '.', snap_pref = 'snapshot'):
		self.dataset = dict.fromkeys(SnapshotContainer.__field_list, [])
		snap = h5py.File('%s/%s_%03d.hdf5' % (snap_dir, snap_pref, snap_no), 'r')
		for key in self.dataset.keys():
			self.dataset[key] = snap[key][()]
			
	def show_field_list(self):
		print(list(self.dataset.keys()))
