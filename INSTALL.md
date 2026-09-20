# Installation instructions for levioSAM2

levioSAM2 supports a variety of methods for installation:

- Pixi (recommended for a reproducible development environment)
- Conda
- Docker
- Singularity
- CMake

## Set up a development environment with Pixi

The checked-in Pixi workspace supports Linux on x86_64 and ARM64, and macOS
on Intel and Apple Silicon. It installs the build tools and HTSlib, then
installs levioSAM2 into the project environment. Linux and Intel macOS use
the conda-forge `sdsl-lite` package. Apple Silicon builds SDSL 2.1.1 from its
pinned source revision because conda-forge does not currently publish an
`osx-arm64` package.

Install a C and C++ compiler first: GCC on Linux or the Xcode Command Line
Tools on macOS. Pixi supplies CMake, Ninja, Git, HTSlib, and the remaining
packaged dependencies.

```shell
pixi run setup
pixi run leviosam2 --help
```

Run the unit tests with:

```shell
pixi run check
```

## Install levioSAM2 via Conda

```shell
conda install -c conda-forge -c bioconda leviosam2
```

Bioconda currently provides native packages for `linux-64`, `linux-aarch64`,
and `osx-64`, but not for Apple Silicon (`osx-arm64`). On Apple Silicon, use
the Pixi instructions above to build and install levioSAM2 natively from source;
the `pixi run setup` task also builds the required SDSL library.

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

We provide a dependency-minimal CTest suite and an extended end-to-end suite.

- Run the unit and lightweight CLI integration tests with `pixi run test`.

- Run the extended pysam and Picard tests with
  `pixi run -e integration test-integration`. Pixi installs the additional
  Python, pysam, and Picard dependencies in a separate environment. The tests
  use temporary output directories and do not modify `testdata`.

- If dependencies are installed manually, run the extended suite from any
  directory with `python /path/to/leviosam-test.py /path/to/leviosam2`.

## ARM64

LevioSAM2 supports the arm64 architecture. Note that the distributed libraries of the dependencies (namely `htslib` and `sdsl-lite`) on Conda or other package managers might not support arm64. Thus, you might need to build the dependent libraries from source. Both `htslib` and `sdsl-lite` can be built under the arm64 architecture.

Notes:

1. `sdsl-lite==v2.1.1` also installs `gtest` and this can result in errors when building levioSAM2 due to duplicated declarations. A work-around is to remove the `gtest` headers installed along with `sdsl-lite`.

2. There can be dynamic library linking errors when executing levioSAM2. This can be solved by `export DYLD_LIBRARY_PATH=/path/to/libsdsl/:/path/to/libhts/:$DYLD_LIBRARY_PATH`.
