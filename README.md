# stmesh

Stmesh is a library implementing the paper "4D Space-Time Delaunay Meshing for Medical Images" by Panagiotis Foteinos and Nikos Chrisochoides. It was devised as part of my Bachelor's Thesis at the RWTH Aachen University.

Through the use of CGAL's efficient implementation of Delaunay triangulations, stmesh is able to run significantly faster than the original implementation. The library is written in C++ and offers the `stmesher` command line tool to mesh 4D geometries.

## Building

This project is designed to be used with devcontainers and in particular works best in VSCode.

To build the project, firstly install the devcontainer extension and open the project in VSCode. You should then be prompted to open the project in a devcontainer. 

Once the devcontainer is running, you can build the project by first configuring the project with CMake:
```bash
cmake --preset unixlike-clang-debug
```

Then build the project with:
```bash
cmake --build --preset unixlike-clang-debug --parallel
```

The `stmesher` binary will be located in the `out/build/unixlike-clang-debug/src/stmesher` directory.

## Running

In order to run `stmesher` on HPC systems, which may sometimes be necessary due to its memory requirement, it is advised to use `apptainer` for containerization. As only a `Dockerfile` is supplied for use with the devcontainer, the docker container must first be built and then converted to a SIF file:
```bash
docker build -t stmesh:latest --build-arg VARIANT=noble --build-arg GCC_VER=14 --build-arg LLVM_VER=18 . --target release
apptainer build stmesh.sif docker-daemon://stmesh:latest
```
This file can then be uploaded to the cluster and used by executing it. This also works in batch files, for more information see the documentation of your cluster.

Note that building directly on the cluster is not possible, as docker is not available there.

## Usage

```
Usage: stmesher [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --version                   Show version information
  --statistics-output TEXT    Write statistics to a file
  --vtk-output-dir TEXT       Directory to write vtk files to
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
  --edt-file TEXT             Read an EDT file
  --no-constant-lfs{false} [1]  Needs: --edt-file
                              Use a constant local feature size, default true
  --hypercube [FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT] Excludes: --use-edt-file-boundary-regions
                              Add a hypercube to the meshing algorithm
  --use-edt-file-boundary-regions Needs: --edt-file Excludes: --hypercube
                              Use boundary regions from the EDT file
  --mixd-output TEXT          Specify the .minf file to write, other MIXD files will be placed alongside it.
```

The `stmesher` tool is used to mesh 4D geometries. By default, it meshes a hypersphere defined by an SDF. This is primarily intended for testing purposes, instead typically it is desirable to specify a binary image as `--edt-file`, which is then used to generate the mesh. It must be supplied in a format that can be read into a `uint8` ITK image.

Boundary regions can be specified in this binary image by through the use of `--use-edt-file-boundary-regions` and using the desired indices as the values in the image at the correct positions. Instead, it is also possible to specify a hypercube in which a specific boundary index will be assigned using the `--hypercube` option.

To get debug information, is is possible to set the `SPDLOG_LEVEL` environment variable to `debug`:
```bash
SPDLOG_LEVEL=debug ./stmesher [...]
```

The parameters `--rho-bar`, `--tau-bar`, `--zeta`, `--b` and `--delta` are used to configure the meshing algorithm. The default values are the ones used in the paper. The `--seed` parameter is used to seed the random number generator, which is used to generate the initial points in the meshing algorithm. `--disable-rule6` can be used to disable the sixth rule of the meshing algorithm, which can cause termination problems.

Additionally, the `--radius-scheme` with `--radius-scheme-arg` parameters can be used to specify the radius scheme used in the meshing algorithm. The default is a `constant` radius scheme, with default infinite radius. Any pentatopes of circumradius larger than this will be split. Also available are the `image` radius scheme reading the radius from an image supplied as the param, the `boundary` radius scheme using the boundary distance as the radius, modified by the expression given in the `--radius-scheme-arg` parameter with `d` as the input variable, and the `lfs` scheme, which instead of the distance to the boundary supplies the distance to the medial axis. The latter two require an EDT file to be supplied.

To calculate the local feature size, one can either choose a constant value, or if an image is used as input, `--no-constant-lfs` can be used to calculate the local feature size from the image via thinning.

For output, the `--vtk-output-dir` parameter can be used to specify a directory to write VTK files to. This can be further configured with the `--vtk-output-name-format`, `--vtk-output-dt`, `--vtk-out-coord-format`, `--vtk-output-vtp-format` and `--vtk-output-blocks` parameters. If the `--vtk-output-vtp-format` options is speicifed, VTP files indicating the boundary indices will be output. Similarly, if `--vtk-out-coord-format` is specified, a binary file in MIXD compatible format containing the coordinates of the mesh will be output.

As VTK output is slow and memory intensive, it may be split into blocks to save memory using the `--vtk-output-blocks` parameter. If it is split into blocks, these may also run in parallel, of course at the cost of some additional memory, by setting `OMP_NUM_THREADS` to the desired number of threads:
```bash
OMP_NUM_THREADS=4 ./stmesher [...]
```

The `--mixd-output` parameter can be used to specify an `.minf` file to write. This file will contain the mesh in a format compatible with XNS. The corresponding `.mxyz`, `.mien` and `.mrng` 
files will be placed alongside the `.minf` file in the same directory. Note that the `.mien` file will not be dualized, however this can be fixed by using the `gendual` tool.

Finally, the `--statistics-output` parameter can be used to write statistics about the final mesh to a CSV file.

The `--output-scale` parameter can be used to scale the output mesh.

## Visualizing simulation results

As the `stmesher` tool is capable of both producing a VTK output and a MIXD output, it is possible to simply add the output of an XNS simulation to the VTU slices. 

For this, firstly the output data must be mapped to the coordinates of the VTU slice vertices. If `--vtk-out-coord-format` was used, the [MeshProjector](https://github.com/JoseAntFer/MeshProjector) tool can be used on it in order to obtain the projected data. For example
```bash
./mesh_projector 4 [stmesh_mxyz_file] [stmesh_mien_file] [cns_output] [stmesh_vtk_out_coord_file] [projected_cns_output] -swap_endianness -fortran_indexing -verbose
```

After this, the data can be added using a simply python script. For instance for a CNS simulation:
```python
import meshio
import numpy as np
from sys import argv

# Usage: python3 add_cns.py <projected_cns_output> <stmesher_vtu> <output_vtu>

dat = np.fromfile(argv[1], dtype='>d').reshape(-1,4)
v = dat[:,:3]
p = dat[:,3]

m = meshio.read(argv[2])
m.point_data["p"] = p
m.point_data["v"] = v

meshio.write(argv[3], m)
```

> ***NOTE:*** Currently there is a bug in meshio cocerning its use with VTU files containing different sized polyhedra. The fix from [this PR](https://github.com/nschloe/meshio/pull/1463) currently needs to be applied manually.

## Developing

With the project open in a devcontainer, you can start developing. In order to get the best experience with VSCode, it is recommended to configure clangd by specifying in `.vscode/settings.json`:
```json
{
    "clangd.arguments": [
        "--compile-commands-dir=out/build/unixlike-clang-debug"
    ]
}
```

### Debugger
Owing to the relatively poor performance in debug mode and potentially rare errors that arise during debugging of the `stmesh` library, it is recommended to set up [`rr`](https://rr-project.org) for reverse debugging support. This can be done by installing `rr` in the devcontainer and then configuring the debugger in `.vscode/launch.json`(example for debugging the tests):
```json
{
    "version": "0.2.0",
    "configurations": [
        {
            "name": "tests (rr/gdb)",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceRoot}/out/build/unixlike-clang-debug/test/tests",
            "args": [],
            "miDebuggerServerAddress": "localhost:50505",
            "stopAtEntry": false,
            "cwd": "${workspaceRoot}/out/build/unixlike-clang-debug/test",
            "environment": [],
            "externalConsole": true,
            "linux": {
                "MIMode": "gdb",
                "setupCommands": [
                    {
                        "description": "Setup to resolve symbols",
                        "text": "set sysroot /",
                        "ignoreFailures": false
                    }
                ]
            }
        }
    ]
}
```
In order to use this, firstly record the execution of the tests with
```bash
rr record out/build/unixlike-clang-debug/test/tests
```
Then start the `rr` server by running the `rr replay -s 50505` command in the terminal. The VSCode debugger can then be started and will connect to the `rr` server. To execute in reverse simply use the debug console and run
```
set exec-direction reverse
```
