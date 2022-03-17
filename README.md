[![Docker](https://img.shields.io/docker/v/naechyun/leviosam2?label=Docker)](https://hub.docker.com/r/naechyun/leviosam2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam2/badges/version.svg)](https://anaconda.org/bioconda/leviosam2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam2/badges/platforms.svg)](https://anaconda.org/bioconda/leviosam2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam2/badges/downloads.svg)](https://anaconda.org/bioconda/leviosam2)

# levioSAM2: Fast and accurate coordinate conversion between assemblies

LevioSAM2 lifts over alignments accurately and efficiently using a chain file.

## Features

- Converting SAM/BAM records from one reference to another using either a VCF file or a chain file
- Converting alignment information, including:
    - Contig
    - Position
    - Alignmant flag
    - CIGAR string
    - Mate read information
    - MD:Z and NM:i tag (optional)

- Multithreading support


## Installation

LevioSAM2 can be installed using:

- [Conda](https://docs.conda.io/en/latest/) (please make sure the version number is specified)

```
conda install -c conda-forge -c bioconda leviosam2
```

- [Docker](https://hub.docker.com/r/milkschen/leviosam2)
```
docker pull naechyun/leviosam2
```

- Built from source using CMake or Make. See [INSTALL.md](INSTALL.md) for details.


## Usage (ChainMap)

LevioSAM2 performs lift-over using a [chain file](http://hgw1.soe.ucsc.edu/goldenPath/help/chain.html) as the lift-over map.

Quick run:
```
leviosam2 index -c a_to_b.chain -p a_to_b -F dest.fai
leviosam2 lift -C a_to_b.clft -a aligned_to_a.bam -p lifted_from_a -O bam
```

The levioSAM2 ChainMap index will be saved to `a_to_b.clft`. The output will be saved to `lifted_from_a.bam`.


## Publication

Taher Mun, Nae-Chyun Chen, Ben Langmead, LevioSAM: Fast lift-over of variant-aware reference alignments, _Bioinformatics_, 2021;, btab396, https://doi.org/10.1093/bioinformatics/btab396
