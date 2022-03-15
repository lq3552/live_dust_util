import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from .SnapshotContainer import SnapshotContainer
from .GrainSizeDistribution import GrainSizeDistribution
from .ExtinctionLaw import ExtinctionLaw


class RadialProfile():

	"""
	Parameters:
		snap: <SnapshotContainer>
		wave: <ndarray[], dtype = float64> wavelength

		op_a: str, filename of grain radii used by optical property tables
		op_gra: str, directory and prefix of optical proprerty table for graphites
		op_sil: str, directory and prefix of optical property table for silicates
		a  : <ndarray[], dtype = float64> centers of grain size bins
		p_c: <ndarray[3], dtype = float64> center to compute radii
		r_s: <float32> lower bound of radius interval
		r_e: <float32> upper bound of radius interval
		nbins: <int> number of bins between r_s and r_e
		lz: <ndarray[3], dtype = float64> direction of angular momentum, default [0, 0, 1]
	"""
	field_list = ["GasDensity",
				  "DustDensity",
				  "SigmaSFR",
				  "BumpStrength",
				  "UVOPSlope",
				  "RV",
				  "CarbonFraction",
				  "ExtinctionCurve"
				 ]

	def __init__(self, snap, wave, op_a = None, op_gra = None, op_sil = None, a = 10**np.linspace(-3,0,16), p_c = [], r_s = None, r_e = None, n_bins = 10, lz = np.array([0.,0.,1.])):
		#may simplize this input list###
		##TODO throw exception if n_bins < 1 or other illegal input
		self.data_radial = dict.fromkeys(RadialProfile.field_list, 0)
		for key in self.data_radial.keys():
			if key != "ExtinctionCurve":
				self.data_radial[key] = np.zeros(n_bins)
			else:
				self.data_radial[key] = np.zeros((n_bins, len(wave)))

		self.rad = np.linspace(r_s, r_e, n_bins)
		if n_bins == 1:
			drad = r_e - r_s
		elif n_bins > 1:
			drad  = self.rad[1] - self.rad[0]
		wave = 10**np.linspace(-1,0,200)
		xyz = snap.dataset["PartType0/Coordinates"]
		xyzd = snap.dataset["PartType3/Coordinates"]
		dr2 = np.sum((xyz - p_c)**2, axis=1)
		drd2 = np.sum((xyzd - p_c)**2, axis=1)

		for i in range(len(self.rad)):
			filt = np.where( (dr2 >= self.rad[i]**2) & (dr2 < (self.rad[i] + drad)**2) )
			filt_d = np.where( (drd2 >= self.rad[i]**2) & (drd2 < (self.rad[i] + drad)**2) )
			gsd = GrainSizeDistribution(snap, a, p_c, self.rad[i], self.rad[i] + drad)
			abu = gsd.compute_abundances()
			ext = ExtinctionLaw(gsd, wave, op_a, op_gra, op_sil)
			try:
				self.data_radial["GasDensity"][i] = np.average(
                                     snap.dataset["PartType0/Density"][filt],
                                     weights = snap.dataset["PartType0/Masses"][filt]
									 )
				self.data_radial["DustDensity"][i] = np.average(
                                     snap.dataset["PartType3/Dust_DustDensity"][filt_d],
                                     weights = snap.dataset["PartType3/Masses"][filt_d]
									 )
				self.data_radial["SigmaSFR"][i] = np.sum(
                                     snap.dataset["PartType0/StarFormationRate"][filt]) / (np.pi* (2 * self.rad[i] * drad + drad**2)) #,
                                     #weights = snap.dataset["PartType0/Masses"][filt]
									 #)
				self.data_radial["BumpStrength"][i] = ext.bump
				self.data_radial["UVOPSlope"][i] = ext.slope_uo
				self.data_radial["RV"][i] = ext.RV
				self.data_radial["CarbonFraction"][i] = (abu['PAH'] + abu['Aliphatic C'])# / abu['Silicate']
				self.data_radial["ExtinctionCurve"][i] = ext.extinction
			except ZeroDivisionError:
				print("RuntimeWarning: find no particles within the annulus [%e, %e)!" 
						% (self.rad[i], self.rad[i] + drad))
