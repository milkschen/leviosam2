# Installation instructions for levioSAM

levioSAM supports a variety of methods for installation:

- Conda (*highly recommended*)
- Docker
- CMake
- Make

## Conda

This is our recommended method of installing levioSAM.

To install, simpy run this command:

```
conda install -c conda-forge -c bioconda leviosam
```

## Docker

You can obtain a Docker image of the latest version from Docker hub:

```
docker pull alshai/leviosam
```

## CMake and Make


Make sure the following prerequisite libraries are installed on your system. 

- [htslib v1.10+](https://github.com/samtools/htslib)
- [sdsl-lite v2.1.1+](https://github.com/simongog/sdsl-lite/)

An easy way to install these dependencies is to use your OS's existing package system:
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
cmake -DHTS_LIB_DIR=<htslib lib directory> \
      -DHTS_INC_DIR=<htslib include dir> \
      -DSDSL_LIB_DIR=<sdsl-lite lib directory> \
      -DSDSL_INC_DIR=<sdsl-lite include dir> \
      ..
```

### Make

Update `LD_LIBRARY_PATH` and `CPLUS_INCLUDE_PATH` paths after installing sdsl-lite and htslib and install with `make`:

```
export LD_LIBRARY_PATH=<path/to/lib>:$LD_LIBRARY_PATH
export C_INCLUDE_PATH=<path/to/include>:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=<path/to/include>:$CPLUS_INCLUDE_PATH
make
```

## Test

We provide an end-to-end test and a set of unit tests for levioSAM.

- The end-to-end test can be run with `python leviosam-test.py`. This test includes running levioSAM on several test files in `testdata`. We also use `picard` to test if the lifted SAM files are valid.

- The unit test can be run `cd build; ctest` if you use cmake to build levioSAM; or `make gtest; cd testdata; ../gtest` if you use make to build levioSAM.
