# [Choice] bionic (18.04), focal (20.04)
ARG VARIANT="focal"
FROM ubuntu:${VARIANT}

# Restate the variant to use it later on in the llvm and cmake installations
ARG VARIANT

# Install necessary packages available from standard repos including python depends
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        software-properties-common wget apt-utils file zip bsdmainutils \
        openssh-client gpg-agent socat rsync \
        make ninja-build git \
        python3 python3-pip python3-venv \
        pre-commit python3-scipy python3-numpy

# User-settable versions:
# This Dockerfile should support gcc-[7, 8, 9, 10, 11] and clang-[10, 11, 12, 13]
ARG GCC_VER="11"
# Add gcc-${GCC_VER}
# temporarily disabled toolchain, not currently required
# RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test 
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        gcc-${GCC_VER} g++-${GCC_VER} gdb

# Set gcc-${GCC_VER} as default gcc
RUN update-alternatives --install /usr/bin/gcc gcc $(which gcc-${GCC_VER}) 100
RUN update-alternatives --install /usr/bin/g++ g++ $(which g++-${GCC_VER}) 100

ARG LLVM_VER="13"
# Add clang-${LLVM_VER}
# temporarily disabled llvm repo, not currently required
#ARG LLVM_URL="http://apt.llvm.org/${VARIANT}/"
#ARG LLVM_PKG="llvm-toolchain-${VARIANT}-${LLVM_VER}"
#RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add - 2>/dev/null && \
#    add-apt-repository -y "deb ${LLVM_URL} ${LLVM_PKG} main"
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

# temporarily use ubuntu version until repository allows noble
# Add current cmake/ccmake, from Kitware
#ARG CMAKE_URL="https://apt.kitware.com/ubuntu/"
#ARG CMAKE_PKG=${VARIANT}
#RUN wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null \
#        | gpg --dearmor - | tee /etc/apt/trusted.gpg.d/kitware.gpg >/dev/null && \
#    apt-add-repository -y "deb ${CMAKE_URL} ${CMAKE_PKG} main" && \
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        cmake cmake-curses-gui

# Install editors
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        neovim emacs nano

# Install optional dependecies
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        doxygen graphviz ccache cppcheck python-is-python3 rr valgrind kcachegrind

# Add required libraries
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        libcgal-dev libvtk9-dev libfftw3-dev patch

# Install mold linker
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        mold

# Add ITK patch to fix compilation error
COPY itk.patch /opt/itk.patch

# Build own version of ITK with patch
RUN git clone --depth 1 --branch v5.2.1 https://github.com/InsightSoftwareConsortium/ITK /opt/itk-src && \
    cd /opt/itk-src && \
    patch -p1 < /opt/itk.patch && \
    mkdir /opt/itk-clang-build && \
    cd /opt/itk-clang-build && \
    CC=clang CXX=clang++ cmake -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF -DITKGroup_Core:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DModule_ITKDistanceMap:BOOL=ON -DITKGroup_IO:BOOL=ON -DModule_ITKCommon:BOOL=ON -DITK_USE_SYSTEM_EIGEN:BOOL=ON /opt/itk-src && \
    cmake --build . -j $(nproc) && \
    mkdir /opt/itk-gcc-build && \
    cd /opt/itk-gcc-build && \
    cmake -DITK_BUILD_DEFAULT_MODULES:BOOL=OFF -DITKGroup_Core:BOOL=OFF -DBUILD_TESTING:BOOL=OFF -DBUILD_EXAMPLES:BOOL=OFF -DModule_ITKDistanceMap:BOOL=ON -DITKGroup_IO:BOOL=ON -DModule_ITKCommon:BOOL=ON -DITK_USE_SYSTEM_EIGEN:BOOL=ON /opt/itk-src && \
    cmake --build . -j $(nproc) && \
    rm /opt/itk.patch

# X11 support
RUN apt-get update -qq && export DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y --no-install-recommends \
        dbus-x11

## Cleanup cached apt data we don't need anymore
RUN apt-get autoremove -y && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Allow the user to set compiler defaults
ARG USE_CLANG
# if --build-arg USE_CLANG=1, set CC to 'clang' or set to null otherwise.
ENV CC=${USE_CLANG:+"clang"}
ENV CXX=${USE_CLANG:+"clang++"}
# if CC is null, set it to 'gcc' (or leave as is otherwise).
ENV CC=${CC:-"gcc"}
ENV CXX=${CXX:-"g++"}

# Include project
#ADD . /workspaces/cpp_starter_project
#WORKDIR /workspaces/cpp_starter_project

CMD ["/bin/bash"]
