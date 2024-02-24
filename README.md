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
cmake --build --preset unixlike-clang-debug --parallel $(nprocs)
```

The `stmesher` binary will be located in the `out/build/unixlike-clang-debug/src/stmesher` directory.

## Running

In order to run `stmesher` on HPC systems, which may sometimes be necessary due to its memory requirement, it is advised to use `apptainer` for containerization. As only a `Dockerfile` is supplied for use with the devcontainer, the docker container must first be built and then converted to a SIF file:
```bash
docker build -t stmesh:latest --build-arg VARIANT=jammy --build-arg GCC_VER=13 --build-arg LLVM_VER=17 .devcontainer
apptainer build stmesh.sif docker-daemon://stmesh:latest
```
This file can then be uploaded to the cluster and used by executing it. This also works in batch files, for more information see the documentation of your cluster.

Note that building directly on the cluster is not possible, as docker is not available there.

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

## Usage

```
Usage: ./stmesher [OPTIONS]

Options:
  -h,--help                   Print this help message and exit
  --version                   Show version information
  --statistics-output TEXT    Write statistics to a file
  --vtk-output-dir TEXT       Directory to write vtk files to
  --vtk-output-name-format TEXT [mesh_{}.vtu]  Needs: --vtk-output-dir
                              Format string for vtk output files
  --vtk-output-dt FLOAT [0.5]  Needs: --vtk-output-dir
                              Time step for vtk output
  --rho-bar FLOAT [20]        Rho bar for meshing algorithm
  --tau-bar FLOAT [0.0013]    Tau bar for meshing algorithm
  --zeta FLOAT [0.5]          Zeta for meshing algorithm
  --b FLOAT [5]               b for meshing algorithm
  --delta FLOAT [5]           delta for meshing algorithm
  --seed UINT                 Seed for random number generation
  --edt-file TEXT             Read an EDT file
  --hypercube [FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT,FLOAT] Excludes: --use-edt-file-boundary-regions
                              Add a hypercube to the meshing algorithm
  --use-edt-file-boundary-regions Needs: --edt-file Excludes: --hypercube
                              Use boundary regions from the EDT file
  --mixd-output TEXT          Specify the .minf file to write, other MIXD files will be placed alongside it.
```
