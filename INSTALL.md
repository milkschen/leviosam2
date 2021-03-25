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
conda install -c bioconda leviosam
```

## Docker

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
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path/to/lib>
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:<path/to/include>
make
```
