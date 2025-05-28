# mesh_analyzer

```
Usage: ./mesh_analyzer [OPTIONS] [SUBCOMMAND]

Options:
  -h,--help                   Print this help message and exit
  --config                    Read an ini file
[Option Group: Version]
  Options:
    --version                   Show version information
[Option Group: Regular]
  Positionals:
    mixd_file_name TEXT REQUIRED
                                The name of the MIXD file to read
  Options:
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
    --data-file TEXT            The name of the data file to read
    --problem-type TEXT:{advection_diffusion,cns,emum,ins,rvtcns,solid,soliduv,viscoelastic} Needs: --data-file
                                The problem type to use for the data file
    --data-entry [TEXT,UINT] ... Needs: --data-file
                                The data entries to construct the problem type from
    --write-config TEXT         Write the config to a file

Subcommands:
  transform                   Apply 4D affine transformation
```

The `mesh_analyzer` tool is used to analyze the quality of 4D meshes in MIXD format and slice them to directly generate time step slices in VTK format, which may optionally include simulation data.

## Input MIXD

As input, a path to an MIXD `.minf` file must be specified. Note that the paths within it are considered relative to the current working directory.

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

#### Data output

Optionally, simulation data may be directly included within the VTK files. For this a data file must be provided via the `--data-file` argument. By default, each entry will be added as a separate point data array named `u_1`,`u_2`,... Instead, the data can be split into more sensible arrays by specifying a known problem type via the `--problem-type` parameter. The available problem types are `advection_diffusion`, `cns`, `emum`, `ins`, `rvtcns`, `solid`, `soliduv` and `viscoelastic`. 

Alternatively, for custom problem types the data entries may be specified directly using the `--data-entry` parameter, which must be specified once for each data entry in order and takes  pairs of the form `name,length`, where `name` is the name of the data entry and `length` is the length of the vector(i.e. 1 for scalars).

#### Transformation subcommand

```
Apply 4D affine transformation
Usage: ./mesh_analyzer transform [OPTIONS]

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


### Statistics output

Finally, the `--statistics-output` parameter can be used to write statistics about the final mesh to a CSV file.

## Debug Information 

To get debug information, is is possible to set the `SPDLOG_LEVEL` environment variable to `debug`:
```bash
SPDLOG_LEVEL=debug ./mesh_analyzer [...]
```

## Configurations

Instead of command line arguments, a configuration file in TOML format may be provided via the `--config` option. In order to generate a configuration file from the given command line arguments, the `--write-config` option may be used. This will write the configuration to the specified file in TOML format to get you started.
