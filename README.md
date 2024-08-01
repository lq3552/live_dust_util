# live\_dust\_util
live_dust_util is a tool to analyze dusty galaxies from `AREPO` simulations with live dust models and `SMUGGLE` stellar formation and feedback model  in `Python`

## Get

You can download it via
```bash
git clone https://github.com/lq3552/live_dust_util.git
```

## Build

live\_dust\_util requires:

 * python(>=3.5)
 * numpy(>=1.12.0)
 * matplotlib(>=3.5.0)
 * [yt](http://yt-project.org)(>=3.4.1)

You can install it by setup.py:
```bash
python -m pip install .
```

or, if you need root access,
```bash
sudo python -m pip install .
```

## Basic usage

To extract extinction curves from GIZMO snapshots:
```python
from live_dust_util import GrainSizeDistribution as GSD
from live_dust_util import ExtinctionLaw as Ext

snap_num = 0
snapshot = Snap(snap_num, "snap_dir")
gsd = GSD(snapshot, a = grain_size_bins, p_c = center_of_galaxy, r_s = radii_lower_bound, r_e = radii_upper_bound)
extinction_law = Ext(gsd, wavelength, 
					 op_a = "gain_size_bin_for_optical_properties.txt",
					 op_gra = "optical_properties_for_graphites.txt",
					 op_sil = "optical_properties_for_silicates.txt"
					)
```
or for snapshots written in parelell, if you want to load a specific sub-snapshot

```python
snap_num = "0".zfill(3)
snap_rank = "1"
snapshot = Snap(snap_num + '.' + snap_rank, "snap_dir")
```
Combination of sub-snapshots is not yet supported.

To analyze dusty galaxy properties or statistics, such as the dust-mass function:
```python
from live_dust_util import Galaxy
gal = Galaxy(snapshot, a = grain_size_bins, p_c = center_of_galaxy, r_s = radii_lower_bound, r_e = radii_upper_bound)
Md = gal.dataset["MassesByType"]["PartType3"]
```
To analyze radial profiles:
```python
import live_dust_util.RadialProfile as RadP
rad_prof = RadP(snap, wavelength,
				a = grain_size_bins,
				op_a = "gain_size_bin_for_optical_properties.txt",
				op_gra = "optical_properties_for_graphites.txt",
				op_sil = "optical_properties_for_silicates.txt"
				p_c = center_of_galaxy,
				r_s = radii_lower_bound,
				r_e = radii_upper_bound,
				n_bins = n_bins
			   )
rad_bins = rad_prof.rad
rad_density = rad_prof.data_radial["GasDensity"]
```

## Examples
For more practical examples, please refer to scripts in `example`
