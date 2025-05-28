# py4dproject

The Python module `py4dproject` is provided in order to simplify evaluating the results of continuous simplex space-time simulations at arbitrary space-time points. It reuses existing facilities within `stmesh` and exposes them in Python via nanobind, and consequently offers good performance characteristics.

## Installation

The module may be extracted from the stmesh docker image and can be found in `/usr/lib/python3.12/site-packages/py4dproject.abi3.so`. Only the Alpine based image contains this module.

## Usage

After importing the module, an instance of the `MeshProjector` class must be created:
```python
from py4dproject import MeshProjector

proj = MeshProjector('artery.minf', 'stokes.out')
```
By default, this will split the data into a dictionary of numpy arrays containing a single column of data, referred to as `u_1`, `u_2`, ...

To instead break the data into more sensible arrays, the `MeshProjector` constructor may be passed a `problem_type` argument, which must be a value of the `ProblemType` enum. The available problem types are `Solid`, `Viscoelastic`, `AdvectionDiffusion`, `INS`, `CNS`, `RVTCNS`, `EMUM` and `SolidUV`:
```python
from py4dproject import MeshProjector, ProblemType

proj = MeshProjector('artery.minf', 'stokes.out', ProblemType.INS)
```

Alternatively, a custom problem type may be specified by passing a list of instances of the DataEntry class to the `MeshProjector` constructor. Each `DataEntry` instance must specify the name of the data entry and its length:
```python
from py4dproject import MeshProjector, DataEntry
proj = MeshProjector(
    'artery.minf',
    'stokes.out',
    [
        DataEntry('velocity', 3),  # 3D velocity vector
        DataEntry('pressure', 1)   # scalar pressure
    ]
)
```

After the `MeshProjector` instance has been created, the `project` method may be used to evaluate the simulation results at a space-time point. The space-time point must be specified as a tuple of `(x, y, z, t)`, where x, y and z are the spatial coordinates and t is the time coordinate:
```python
result = proj.project(coordinates)
```
Note that the data type of the coordinates must be `np.float64` with a dimension of `(n, 4)`.

## Compatibility

`py4dproject` is built with the stable Python3 ABI. Due to limitations of the binder(nanobind), it is however only comatible with Python 3.12 and later.

As the `py4dproject.abi3.so` statically links all required dependencies, including the C++ standard library and `musl-libc`(which is why it can only be built on musl based systems), it can be used on any Linux system regardless of the distribution or version.
