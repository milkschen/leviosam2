[![Docker](https://img.shields.io/docker/v/naechyun/leviosam2?label=Docker)](https://hub.docker.com/r/naechyun/leviosam2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam2/badges/version.svg)](https://anaconda.org/bioconda/leviosam2)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam2/badges/downloads.svg)](https://anaconda.org/bioconda/leviosam2)

# levioSAM2: Fast and accurate coordinate conversion between assemblies

LevioSAM2 lifts over alignments accurately and efficiently using a chain file.

## Features

- Converting aligned short and long reads records (in SAM/BAM format) from one reference to another
- The following alignment features are updated during coversion:
    - Reference name (`RNAME`), position (`POS`), alignmant flag (`FLAG`), and CIGAR alignment string (`CIGAR`)
    - Mate read information (`RNEXT`, `PNEXT`, `TLEN`)
    - (optional) Alignment tags (`MD:Z`, `NM:i`)
- Multithreading support
- Toolkit for "selective" pipelines which consider major changes between the source and target references
- (beta) Converting intervals (in BED format) from one reference to another

## Installation

LevioSAM2 can be installed using:

- [Conda](https://docs.conda.io/en/latest/)

```
conda install -c conda-forge -c bioconda leviosam2
```

- [Docker](https://hub.docker.com/r/milkschen/leviosam2)
```
docker pull naechyun/leviosam2:latest
```

- [Singularity](https://hub.docker.com/r/milkschen/leviosam2)
```
singularity pull docker://naechyun/leviosam2:latest
```

- Built from source using CMake. See [INSTALL.md](INSTALL.md) for details.


## Usage

LevioSAM2 performs lift-over using a [chain file](http://hgw1.soe.ucsc.edu/goldenPath/help/chain.html) as the lift-over map.

Quick run:
```
leviosam2 index -c source_to_target.chain -p source_to_target -F target.fai
leviosam2 lift -C source_to_target.clft -a aligned_to_source.bam -p lifted_from_source -O bam
```

The levioSAM2 ChainMap index will be saved to `source_to_target.clft`. The output will be saved to `lifted_from_source.bam`.

See [this README](https://github.com/milkschen/leviosam2/blob/main/pipeline/README.md) to learn more about running the full levioSAM2
pipeline with selective re-mapping.


## Publication

Taher Mun, Nae-Chyun Chen, Ben Langmead, LevioSAM: Fast lift-over of variant-aware reference alignments, _Bioinformatics_, 2021;, btab396, https://doi.org/10.1093/bioinformatics/btab396
