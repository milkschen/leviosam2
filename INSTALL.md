# Installation instructions for levioSAM2

levioSAM2 supports a variety of methods for installation:

- Conda
- Docker
- Singularity
- CMake

## Conda

```shell
conda install -c conda-forge -c bioconda leviosam2
```

## Docker

You can obtain a Docker image of the latest version from Docker hub:

```shell
docker pull naechyun/leviosam2:latest
```

## Singularity

```shell
singularity pull docker://naechyun/leviosam2:latest
```

## CMake

### Prerequisites

Make sure the following prerequisite libraries are installed on your system. 

- [htslib v1.12+ (tested up to v1.18)](https://github.com/samtools/htslib)
- [sdsl-lite v2.1.1+](https://github.com/simongog/sdsl-lite/)

Both libraries are available through coda:

```shell
conda install -c conda-forge sdsl-lite
conda install -c bioconda htslib
```

Another easy way to install these dependencies is to use your OS's existing package system:

```shell
apt-get install libhts-dev libsdsl-dev # Debian/Ubuntu
brew tap brewsci/bio; brew install htslib sdsl-lite # MacOS
```

If using RedHat or Fedora, then you must install sdsl-lite manually. But you can install htslib through yum:

```shell
yum install htslib
```

Or you can choose to install them manually by following the install instructions on their respective pages.

### CMake

Once the prerequisite packages are installed and the locations of their installations are known, specify their locations
to CMake by running the following commands:

```shell
mkdir build
cd build
cmake ..
make
# make install
# or 
# make install DESTDIR=/path/to/install
```

If you installed the dependencies manually, you might have to modify the `cmake` command to specify their library and
include directory locations like so:

```shell
cmake -D CMAKE_LIBRARY_PATH="/path/to/libsdsl/;/path/to/libhts/" \
      -D CMAKE_INCLUDE_PATH="/path/to/libsdsl/include/;/path/to/libhts/include/" ..
```

## Test

We provide an end-to-end test and a set of unit tests for levioSAM.

- The end-to-end test can be run with `python leviosam-test.py`. This test includes running levioSAM on several test files in `testdata`. We also use `picard` to test if the lifted SAM files are valid.

- The unit test can be run `cd build; ctest` if you use cmake to build levioSAM; or `make gtest; cd testdata; ../gtest` if you use make to build levioSAM2.

## ARM64

LevioSAM2 supports the arm64 architecture. Note that the distributed libraries of the dependencies (namely `htslib` and `sdsl-lite`) on Conda or other package managers might not support arm64. Thus, you might need to build the dependent libraries from source. Both `htslib` and `sdsl-lite` can be built under the arm64 architecture.

Notes:

1. `sdsl-lite==v2.1.1` also installs `gtest` and this can result in errors when building levioSAM2 due to duplicated declarations. A work-around is to remove the `gtest` headers installed along with `sdsl-lite`.

2. There can be dynamic library linking errors when executing levioSAM2. This can be solved by `export DYLD_LIBRARY_PATH=/path/to/libsdsl/:/path/to/libhts/:$DYLD_LIBRARY_PATH`.
