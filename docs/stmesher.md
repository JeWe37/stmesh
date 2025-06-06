# stmesher

```
Usage: ./stmesher [OPTIONS] [SUBCOMMAND]

Options:
  -h,--help                   Print this help message and exit
  --version                   Show version information
  --statistics-output TEXT    Write statistics to a file
  --vtk-output-dir TEXT:DIR   Directory to write vtk files to
  --vtk-output-name-format TEXT [mesh_{}.vtu]  Needs: --vtk-output-dir
                              Format string for vtk output files
  --vtk-output-dt FLOAT [0.5]  Needs: --vtk-output-dir
                              Time step for vtk output
  --vtk-out-coord-format TEXT Needs: --vtk-output-dir
                              Format string for the coordinates file name
  --vtk-output-vtp-format TEXT Needs: --vtk-output-dir
                              Format string for the vtp output files
  --vtk-output-blocks UINT [1]  Needs: --vtk-output-dir
                              Number of blocks for vtk output
  --output-scale FLOAT [1]    Scale factor for output
  --rho-bar FLOAT [20]        Rho bar for meshing algorithm
  --tau-bar FLOAT [0.0013]    Tau bar for meshing algorithm
  --zeta FLOAT [0.5]          Zeta for meshing algorithm
  --b FLOAT [5]               b for meshing algorithm
  --delta FLOAT [5]           delta for meshing algorithm
  --radius-scheme TEXT [constant] 
                              Radius scheme for meshing algorithm
  --radius-scheme-arg TEXT [inf] 
                              Argument for the radius scheme
  --disable-rule6             Disable picking region
  --seed UINT                 Seed for random number generation
  --hypercube [FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT] ... Excludes: --use-edt-file-boundary-regions
                              Add a hypercube to the meshing algorithm
  --mixd-output TEXT          Specify the .minf file to write, other MIXD files will be placed alongside it.
  --ideal-mixd-positions Needs: --mixd-output
                              Write ideal positions to the MIXD file
  --compute-dual Needs: --mixd-output
                              Include the dual in the mrng file
  --write-neim Needs: --mixd-output
                              Write the node element index mapping
  --write-config TEXT         Write the config to a file
  --config                    Read an ini file

Subcommands:
  edt_geometry                Use an EDT file for geometry
  hypercube_geometry          Use hypercube as geometry
  hypersphere_geometry        Use hypersphere as geometry
  cylinder_geometry           Use time-extruded cylinder as geometry
  transform                   Apply 4D affine transformation
```

The `stmesher` tool is used to mesh 4D geometries. By default, it meshes a hypersphere defined by an SDF. This is primarily intended for testing purposes, instead typically a geometry may be specified via one of the subcommand. 

## Geometry specification

### SDF Geometries

Analytical geometries can be defined via SDFs (Signed Distance Functions). These are faster and accurate, but cannot be used for general geometries.

#### Hypercube

```
Use hypercube as geometry
Usage: ./stmesher hypercube_geometry [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --hypercube-sizes [FLOAT,FLOAT,FLOAT,FLOAT] REQUIRED
                              Sizes of the hypercube
```

A hypercube of dimensions given by the `--hypercube-sizes` parameter is used as the geometry.

#### Hypersphere

```
Use hypersphere as geometry
Usage: ./stmesher hypersphere_geometry [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --hypersphere-center [FLOAT,FLOAT,FLOAT,FLOAT]
                              Center of the hypersphere
  --hypersphere-radius FLOAT  Radius of the hypersphere
```

A hypersphere of radius given by the `--hypersphere-radius` parameter and center given by the `--hypersphere-center` parameter is used as the geometry. Note that typically the center should be located at `t=r`.

#### Cylinder

```
Use time-extruded cylinder as geometry
Usage: ./stmesher cylinder_geometry [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --cylinder-radius FLOAT REQUIRED
                              Radius of the cylinder
  --cylinder-height FLOAT REQUIRED
                              Height of the cylinder
  --time-extrusion FLOAT REQUIRED
                              Time extrusion of the cylinder
```

A cylinder of radius given by the `--cylinder-radius` parameter, height given by the `--cylinder-height` parameter and time extrusion given by the `--time-extrusion` parameter is used as the geometry. The cylinder is extruded in time, so that it has a height in that direction of `time-extrusion`.

### Image geometry

```
Use an EDT file for geometry
Usage: ./stmesher edt_geometry [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --edt-file TEXT REQUIRED    Read an EDT file
  --no-constant-lfs{false}    Use a constant local feature size, default true
  --use-edt-file-boundary-regions Excludes: --hypercube
                              Use boundary regions from the EDT file
```

The geometry is specified as a marching-cubes surface from a binary image. It must be supplied in a format that can be read into a `uint8` ITK image to the requried `--edt-file` parameter. 

Using this geometry, instead of using hypercubes, the `--use-edt-file-boundary-regions` option can be used to specify boundary regions, which will be assigned indices depending on the value present in the image file.

If the `--no-constant-lfs` option is specified, the local feature size will be calculated from the image via thinning. Otherwise, a constant local feature size will be used, which can be set using the `--delta` argument as usual.

## Meshing algorithm configuration

The parameters `--rho-bar`, `--tau-bar`, `--zeta`, `--b` and `--delta` are used to configure the meshing algorithm. The default values are the ones used in the paper. The `--seed` parameter is used to seed the random number generator, which is used to generate the initial points in the meshing algorithm. `--disable-rule6` can be used to disable the sixth rule of the meshing algorithm, which can cause termination problems.

### Radius schemes

Additionally, the `--radius-scheme` with `--radius-scheme-arg` parameters can be used to specify the radius scheme used in the meshing algorithm. The default is a `constant` radius scheme, with default infinite radius. Any pentatopes of circumradius larger than this will be split. Also available are the `image` radius scheme reading the radius from an image supplied as the param, the `boundary` radius scheme using the boundary distance as the radius, modified by the expression given in the `--radius-scheme-arg` parameter with `d` as the input variable, and the `lfs` scheme, which instead of the distance to the boundary supplies the distance to the medial axis. The latter two require an EDT file to be supplied.

### LFS schemes

To calculate the local feature size, one can either choose a constant value, or if an image is used as input, `--no-constant-lfs` can be used to calculate the local feature size from the image via thinning.

## Output formats

Boundary regions, which can be used to set boundary conditions in simulations, can be specified as hypercubes in which a specific boundary index will be assigned using the `--hypercube` option.

Alternatively, if using an EDT image as input, the `--use-edt-file-boundary-regions` option may be used which will assign the indices depending on the value present in the image file.

The `--output-scale` parameter can be used to scale the output mesh.

### VTK output

For output, the `--vtk-output-dir` parameter can be used to specify a directory to write VTK files to. This can be further configured with the `--vtk-output-name-format`, `--vtk-output-dt`, `--vtk-out-coord-format`, `--vtk-output-vtp-format` and `--vtk-output-blocks` parameters. If the `--vtk-output-vtp-format` options is speicifed, VTP files indicating the boundary indices will be output. Similarly, if `--vtk-out-coord-format` is specified, a binary file in MIXD compatible format containing the coordinates of the mesh will be output.

As VTK output is slow and memory intensive, it may be split into blocks to save memory using the `--vtk-output-blocks` parameter. If it is split into blocks, these may also run in parallel, of course at the cost of some additional memory, by setting `OMP_NUM_THREADS` to the desired number of threads:
```bash
OMP_NUM_THREADS=4 ./stmesher [...]
```

#### Transformation subcommand

```
Apply 4D affine transformation
Usage: ./stmesher transform [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  -t,--translate [FLOAT,FLOAT,FLOAT,FLOAT] x 4 Excludes: --matrix
                              4D translation vector (x,y,z,w)
  --rotate-xy FLOAT Excludes: --matrix
                              Rotation angle in radians for XY plane
  --rotate-xz FLOAT Excludes: --matrix
                              Rotation angle in radians for XZ plane
  --rotate-xw FLOAT Excludes: --matrix
                              Rotation angle in radians for XW plane
  --rotate-yz FLOAT Excludes: --matrix
                              Rotation angle in radians for YZ plane
  --rotate-yw FLOAT Excludes: --matrix
                              Rotation angle in radians for YW plane
  --rotate-zw FLOAT Excludes: --matrix
                              Rotation angle in radians for ZW plane
  -s,--scale [FLOAT,FLOAT,FLOAT,FLOAT] x 4 Excludes: --matrix
                              Scale factors (x,y,z,w)
  --matrix FLOAT x 20 Excludes: --translate --rotate-xy --rotate-xz --rotate-xw --rotate-yz --rotate-yw --rotate-zw --scale
                              Custom transform matrix (20 elements, row-major)
```

The `transform` subcommand can be used to apply a 4D affine transformation to the mesh. This is useful for transforming the mesh into a different coordinate system before outputting it. This for instance allows slicing to be performed along different axes.

The transformation may either be specified as a combination of a non-uniform scaling, translation and rotation via the `--translate`, `--rotate-xy`, `--rotate-xz`, `--rotate-xw`, `--rotate-yz`, `--rotate-yw`, `--rotate-zw` and `--scale` parameters, or as a custom transformation matrix via the `--matrix` parameter. The latter is specified in row-major order, with 20 elements for a 4D affine transformation.


### MIXD output

The `--mixd-output` parameter can be used to specify an `.minf` file to write. This file will contain the mesh in a format compatible with XNS. The corresponding `.mxyz`, `.mien` and `.mrng` files will be placed alongside the `.minf` file in the same directory.

For easier partitioning, the `--compute-dual` option may be passed, which will cause a dualized `.mien` file to be generated. Additionally, a `.neim` file may be written using the `--write-neim` option, which contains the node-element index mapping.

Finally, the `--ideal-mixd-positions` option may be used to write the ideal positions of the nodes to the MIXD file. These can be used to move boundary vertices closer to their ideal positions on the geometry boundary in external post-processing.

### Statistics output

Finally, the `--statistics-output` parameter can be used to write statistics about the final mesh to a CSV file.

## Debug Information 

To get debug information, is is possible to set the `SPDLOG_LEVEL` environment variable to `debug`:
```bash
SPDLOG_LEVEL=debug ./stmesher [...]
```

## Configurations

Instead of command line arguments, a configuration file in TOML format may be provided via the `--config` option. In order to generate a configuration file from the given command line arguments, the `--write-config` option may be used. This will write the configuration to the specified file in TOML format to get you started.
