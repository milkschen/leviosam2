# Installation instructions for levioSAM2

levioSAM2 supports a variety of methods for installation:

- Conda
- Docker
- Singularity
- CMake

## Install levioSAM2 via Conda

```shell
conda install -c conda-forge -c bioconda leviosam2
```

## Use a levioSAM2 Docker image

You can obtain a Docker image of the latest version from Docker hub:

```shell
docker pull naechyun/leviosam2:latest
```

## Use a levioSAM2 Singularity image

```shell
singularity pull docker://naechyun/leviosam2:latest
```

## Intall levioSAM2 from scratch using CMake

### Dependencies

Make sure the following prerequisite libraries are installed on your system.

- [htslib v1.12+ (tested up to v1.18)](https://github.com/samtools/htslib)
- [sdsl-lite v2.1.1+](https://github.com/simongog/sdsl-lite/)

The dependent libraries can be installed through the following package manager options, or built from scratch:

```shell
# Conda
conda install -c conda-forge sdsl-lite
conda install -c bioconda htslib

# Debian/Ubuntu
apt-get install libhts-dev libsdsl-dev 

# MacOS
brew tap brewsci/bio
brew install htslib sdsl-lite

# RedHat or Fedora
yum install htslib
# sdsl-lite needs to be installed manually
```

### CMake

Command:

```shell
mkdir build
cd build
cmake ..
make
# make install
# or 
# make install DESTDIR=/path/to/install
```

If you installed the dependencies manually, you might need to provide the path
of the dependent libraries to `cmake` by using the following command:

```shell
cmake -D CMAKE_LIBRARY_PATH="/path/to/libsdsl/;/path/to/libhts/" \
      -D CMAKE_INCLUDE_PATH="/path/to/libsdsl/include/;/path/to/libhts/include/" ..
```

## Testing

We provide an end-to-end test and a set of unit tests for levioSAM.

- The end-to-end test can be run with `python leviosam-test.py`. This test includes running levioSAM on several test files in `testdata`. We also use `picard` to validate the lifted results.

- The unit test can be run with the command `cd build; ctest --verbose`.

## ARM64

LevioSAM2 supports the arm64 architecture. Note that the distributed libraries of the dependencies (namely `htslib` and `sdsl-lite`) on Conda or other package managers might not support arm64. Thus, you might need to build the dependent libraries from source. Both `htslib` and `sdsl-lite` can be built under the arm64 architecture.

Notes:

1. `sdsl-lite==v2.1.1` also installs `gtest` and this can result in errors when building levioSAM2 due to duplicated declarations. A work-around is to remove the `gtest` headers installed along with `sdsl-lite`.

2. There can be dynamic library linking errors when executing levioSAM2. This can be solved by `export DYLD_LIBRARY_PATH=/path/to/libsdsl/:/path/to/libhts/:$DYLD_LIBRARY_PATH`.
