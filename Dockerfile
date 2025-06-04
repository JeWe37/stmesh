# [Choice] bionic (18.04), focal (20.04)
ARG VARIANT="noble"
FROM ubuntu:${VARIANT} AS build

# Restate the variant to use it later on in the llvm and cmake installations
ARG VARIANT

# Install base build packages
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        software-properties-common wget apt-utils file zip bsdmainutils \
        make ninja-build git \
        python3 python3-pip python3-venv \
        python3-scipy python3-numpy python-is-python3

# Install GCC toolchain
ARG GCC_VER="14"
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        gcc-${GCC_VER} g++-${GCC_VER}

RUN update-alternatives --install /usr/bin/gcc gcc $(which gcc-${GCC_VER}) 100
RUN update-alternatives --install /usr/bin/g++ g++ $(which g++-${GCC_VER}) 100

# Add CMAKE packages
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        cmake cmake-curses-gui

# Install mold linker
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        mold

# Add required libraries
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        libcgal-dev libvtk9-dev libfftw3-dev patch

# Add ITK patch to fix compilation error
COPY docker/itk.patch /opt/itk.patch

ARG CORES

# Build own version of ITK with patch
RUN git clone --depth 1 --branch v5.2.1 https://github.com/InsightSoftwareConsortium/ITK /opt/itk-src && \
    cd /opt/itk-src && \
    patch -p1 < /opt/itk.patch && \
    mkdir /opt/itk-gcc-build && \
    cd /opt/itk-gcc-build && \
    cmake -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF -DITKGroup_Core:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DModule_ITKDistanceMap:BOOL=ON -DITKGroup_IO:BOOL=ON -DModule_ITKCommon:BOOL=ON -DITK_USE_SYSTEM_EIGEN:BOOL=ON /opt/itk-src && \
    cmake --build . -j ${CORES:-$(nproc)}

## Cleanup cached apt data we don't need anymore
RUN apt-get autoremove -y && apt-get clean

FROM build AS release

# Build the project
WORKDIR /stmesh
RUN --mount=type=bind,target=.,src=.,rw cmake --preset unixlike-gcc-release && \
                                        cmake --build --preset unixlike-gcc-release -j ${CORES:-$(nproc)} && \
                                        cmake --install out/build/unixlike-gcc-release/ --prefix /usr

# Cleanup build files
VOLUME /out
WORKDIR /out
RUN rm /opt/stmesh -rf

ENTRYPOINT ["/usr/bin/stmesher"]

FROM build AS develop

# Install clang toolchain
ARG LLVM_VER="18"
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        clang-${LLVM_VER} lldb-${LLVM_VER} lld-${LLVM_VER} clangd-${LLVM_VER} \
        llvm-${LLVM_VER}-dev libclang-${LLVM_VER}-dev clang-tidy-${LLVM_VER} \
        clang-format-${LLVM_VER} libclang-rt-${LLVM_VER}-dev libomp-${LLVM_VER}-dev

# Set the default clang-tidy, clang-format and clangd, so CMake/VSCode can find it
RUN update-alternatives --install /usr/bin/clang-tidy clang-tidy $(which clang-tidy-${LLVM_VER}) 1
RUN update-alternatives --install /usr/bin/clang-format clang-format $(which clang-format-${LLVM_VER}) 1
RUN update-alternatives --install /usr/bin/clangd clangd $(which clangd-${LLVM_VER}) 1

# Set clang-${LLVM_VER} as default clang
RUN update-alternatives --install /usr/bin/clang clang $(which clang-${LLVM_VER}) 100
RUN update-alternatives --install /usr/bin/clang++ clang++ $(which clang++-${LLVM_VER}) 100

# Install dev packages
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        openssh-client gpg-agent socat rsync \
        pre-commit gdb

# Install editors
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        neovim emacs nano

# Install optional dependecies
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        doxygen graphviz ccache cppcheck rr valgrind kcachegrind

# Build clang version of ITK
RUN mkdir /opt/itk-clang-build && \
    cd /opt/itk-clang-build && \
    CC=clang CXX=clang++ cmake -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF -DITKGroup_Core:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DModule_ITKDistanceMap:BOOL=ON -DITKGroup_IO:BOOL=ON -DModule_ITKCommon:BOOL=ON -DITK_USE_SYSTEM_EIGEN:BOOL=ON /opt/itk-src && \
    cmake --build . -j ${CORES:-$(nproc)}

# X11 support
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        dbus-x11

## Cleanup cached apt data we don't need anymore
RUN apt-get autoremove -y && apt-get clean

# Allow the user to set compiler defaults
ARG USE_CLANG
# if --build-arg USE_CLANG=1, set CC to 'clang' or set to null otherwise.
ENV CC=${USE_CLANG:+"clang"}
ENV CXX=${USE_CLANG:+"clang++"}
# if CC is null, set it to 'gcc' (or leave as is otherwise).
ENV CC=${CC:-"gcc"}
ENV CXX=${CXX:-"g++"}

CMD ["/bin/bash"]
