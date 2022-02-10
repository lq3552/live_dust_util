from distutils.core import setup
from distutils import util
from live_dust_util import __version__ as version

path_utils = util.convert_path('live_dust_util/utils')

setup(name='live_dust_util',
		version=version,
		description = 'A tool to analyze dust properties in galaxies from Arepo_dust snapshots',
		author = 'Qi Li',
		author_email=['lq3552@gmail.com'], 
		url='https://bitbucket.org/lq3552/live_dust_util',
		requires=['numpy(>=1.12.0)', 'matplotlib(>=2.0.0)','scipy(>=1.6.2)','h5py(>=2.10.0)'],
		package_dir={'live_dust_util': 'live_dust_util', 'live_dust_util.utils': path_utils}, # the present directory maps to src 
		packages = ['live_dust_util', 'live_dust_util.utils']
	 )
