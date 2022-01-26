from dataclasses import dataclass
from .Galaxy import Galaxy

@dataclass
class GalaxyCatalog(list):
	"""
	This class is currently empty in essence and only used by myself to 
	play with inheritance and decorators.
	Just wait for our up-coming cosmological simulations and see!
	"""

	def add(self, galaxy):
		self.append(galaxy)

	def remove(self):
		pass

	def print_info(self):
		pass
