# Installation instructions for levioSAM2

levioSAM2 supports a variety of methods for installation:

- Conda
- Docker
- Singularity
- CMake

## Conda

```
conda install -c conda-forge -c bioconda leviosam2
```

## Docker

You can obtain a Docker image of the latest version from Docker hub:

```
docker pull naechyun/leviosam2:v0.1.0
```

## Singularity

```
singularity pull docker://naechyun/leviosam2:v0.1.0
```

## CMake

### Prerequisites

Make sure the following prerequisite libraries are installed on your system. 

- [htslib v1.10+](https://github.com/samtools/htslib)
- [sdsl-lite v2.1.1+](https://github.com/simongog/sdsl-lite/)

Both libraries are available through coda:
```
conda install -c conda-forge sdsl-lite
conda install -c bioconda htslib
```

Another easy way to install these dependencies is to use your OS's existing package system:
```
apt-get install libhts-dev libsdsl-dev # Debian/Ubuntu
brew tap brewsci/bio; brew install htslib sdsl-lite # MacOS
```

If using RedHat or Fedora, then you must install sdsl-lite manually. But you can install htslib through yum:
```
yum install htslib
```

Or you can choose to install them manually by following the install instructions on their respective pages.

### CMake 

Once the prerequisite packages are installed and the locations of their installations are known, specify their locations
to CMake by running the following commands:

```
mkdir build
cd build
cmake ..
make
```

If you installed the dependencies manually, you might have to modify the `cmake` command to specify their library and
include directory locations like so:

```
cmake -D CMAKE_LIBRARY_PATH="/path/to/libsdsl/;/path/to/libhts/" \
      -D CMAKE_INCLUDE_PATH="/path/to/include/" ..
```


## Test

We provide an end-to-end test and a set of unit tests for levioSAM.

- The end-to-end test can be run with `python leviosam-test.py`. This test includes running levioSAM on several test files in `testdata`. We also use `picard` to test if the lifted SAM files are valid.

- The unit test can be run `cd build; ctest` if you use cmake to build levioSAM; or `make gtest; cd testdata; ../gtest` if you use make to build levioSAM.
