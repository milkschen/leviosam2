# levioSAM lifts alignments to the reference genome

Use a VCF file containing alternative haplotype information to lift SAM/BAM alignments
from that haplotype to the reference sequence

## Building

First, clone the repository:

```
git clone https://github.com/alshai/levioSAM
```

### Dependencies

- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [htslib](https://github.com/samtools/htslib)

You can install these dependencies using conda:

```
conda install -c conda-forge sdsl-lite
conda install -c bioconda htslib
```

Or homebrew (for Macs)
```
brew install htslib
brew install sdsl-lite
```

[CMake](https://cmake.org) or `make` can be used to build `levioSAM`.


### Using CMake

Our `CMakeLists.txt` file expectes to use `PkgConfig` to load the libraries; you may need to add the `pkgconfig`
subdirectories of the `sdsl-lite` and `libhts` libraries to the `CMAKE_PREFIX_PATH` environment variable yourself.

```
mkdir build
cd build
cmake ..
make
```

### Using make

Update `LD_LIBRARY_PATH` and `CPLUS_INCLUDE_PATH` paths after installing sdsl-lite and htslib and install with `make`:

```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:<path/to/lib>
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:<path/to/include>
make
```


## Usage (command line)

**DISCLAIMER**: normalize and left-align your VCF files before lifting alignments!

- We suggest using `bcftools norm` to do this.

To serialize lift-over information between a reference sequence and given alternative sequence described in a VCF file:

```
$ ./levioSAM serialize -v <vcf> -s <sample_name> -p <output prefix>
```

The levioSAM file saved to `<output prefix>.lft`.

To lift over coordinates given from a SAM/BAM file using a serialized `.lft` file:

```
$ ./levioSAM lift -a <sam> -l <lft> -p <output prefix>
```

The lifted coordinates will be saved to `<output prefix>.sam`.

To lift over coordinates without serializing (note: this will be slower):

```
$ ./levioSAM lift -a <sam> -v <vcf> -s <sample_name> -p <output prefix>
```

To read from stdin, use `-a -` or exclude `-a`.


## Usage (C++)

Use the `LiftMap` class to generate lift-over information between a reference genome and an alternative genotype.
A `VCF` file containing the `FMT/GT` field for a specified sample must be provided.

```
#include <cstdio>
#include <htslib/vcf.h>
#include "levioSAM.hpp"

int main() {
    const char* fname = "data/dna.vcf";
    vcfFile* fp = bcf_open(fname, "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    lift::LiftMap l(fp, hdr, "<sample_name>");
}
```

To query the equivalent reference position for a given haplotype position:

```
    l.alt_to_ref("contig_name", 8)); // give the contig name and the position on the contig
```

To serialize lift-over information to a file:

```
    std::ofstream o("alt_to_ref.lft");
    l.serialize(o);
    o.close();
```

To load from a serialized file

```
    std::ifstream in("alt_to_ref.lft");
    lift::Lift l2(in);
    in.close();
```


## Features Currently Support:

- Serialized lift-over information VCF file w/ FMT/GT field for a specified sample.
- Convert SAM/BAM records from haplotype to reference.
- Multithreading support.
- Recalculate read-pair information.
- Recalculate `MD:z` and `NM:i` tags with `-m` (optionally).

## TODO

- Recalculate MAPQ
- Support chain format
