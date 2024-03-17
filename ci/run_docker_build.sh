#!/bin/sh
docker build -t stmesh-builder --build-arg='VARIANT=jammy' --build-arg='GCC_VER=13' --build-arg='LLVM_VER=17' .devcontainer
case $1 in
    build)
        docker run --user $(id -u ${USER}):$(id -g ${USER}) -v $PWD:/repo -w /repo stmesh-builder /repo/ci/run_build.sh $2
        ;;
    test)
        docker run --user $(id -u ${USER}):$(id -g ${USER}) -v $PWD:/repo -w /repo stmesh-builder ctest --preset test-$2
        ;;
esac
