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
python setup.py install
```

or, if you need root access,
```bash
sudo python setup.py install
```

## Basic usage

To extract extinction curves from GIZMO snapshots:
```python
from live_dust_util import SnapshotContainer as Snap
from live_dust_util import ExtinctionLaw
snapshot = Snap('snapshot_000.hdf5')
extinction_law = ExtinctionLaw(snapshot)
```

To analyze dusty galaxy properties or statistics, such as the dust-mass function:
```python
from live_dust_util import Galaxy
from live_dust_util import SnapshotContainer as Snap
snapshot = Snap('snapshot_000.hdf5')
galaxy = Galaxy(snapshot)
```

## Advanced atrributes

```python
# under construction
```
