
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam/badges/installer/conda.svg)](https://conda.anaconda.org/bioconda)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam/badges/version.svg)](https://anaconda.org/bioconda/leviosam)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam/badges/platforms.svg)](https://anaconda.org/bioconda/leviosam)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/leviosam/badges/downloads.svg)](https://anaconda.org/bioconda/leviosam)

# levioSAM: accurate and efficient lift-over of alignments

LevioSAM lifts over alignments accurately and efficiently using a VCF or a chain file.

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

LevioSAM can be installed using:

- [Conda](https://docs.conda.io/en/latest/) (please make sure the version number is specified)

```
conda install -c conda-forge -c bioconda leviosam=0.5.1
```

- [Docker](https://hub.docker.com/r/alshai/leviosam)
```
docker pull alshai/leviosam
```

- Built from source using CMake or Make. See [INSTALL.md](INSTALL.md) for details.


## Usage (VcfMap)

LevioSAM can perform lift-over using a VCF file as the lift-over map. Please visit the [VcfMap Usage Page](https://github.com/alshai/levioSAM/wiki/Lift-over-using-a-VCF-map) for detailed instructions.

Quick run:
```
leviosam index -v <vcf> -s <sample_name> -p <output_prefix>
leviosam lift -a <sam> -l <lft> -p <output_prefix> -O bam
```
The levioSAM VcfMap index will be saved to `<output_prefix>.lft`. The output will be saved to `<output_prefix>.bam`.

### Example (with pre-built 1000-Genomes-Project indexes)

We provide the Bowtie 2 (which are compatible with Bowtie, too!) and levioSAM indexes for the major-allele references based on the 1000 Genomes Project. 
Please navigate to the [bowtie-majref](https://github.com/BenLangmead/bowtie-majref) repo for more comprehensive description and more resources.

We provide more detailed instructions of how to use levioSAM in common variant-aware reference pipelines (major-allele reference and personalized reference) in the [levioSAM wiki](https://github.com/alshai/levioSAM/wiki/Alignment-with-variant-aware-reference-genomes). 


## Usage (ChainMap)

LevioSAM can also perform lift-over using a chain file as the lift-over map since v0.5.0. Using ChainMap enables assembly-to-assembly lift-over and is usually faster. Please visit the [ChainMap Usage Page](https://github.com/alshai/levioSAM/wiki/Lift-over-using-a-chain-map) for detailed instructions.

Quick run:
```
leviosam index -c a_to_b.chain -p a_to_b -F dest.fai
leviosam lift -C a_to_b.clft -a aligned_to_a.bam -p lifted_from_a -O bam
```

The levioSAM ChainMap index will be saved to `a_to_b.clft`. The output will be saved to `lifted_from_a.bam`.


## Publication

Taher Mun, Nae-Chyun Chen, Ben Langmead, LevioSAM: Fast lift-over of variant-aware reference alignments, _Bioinformatics_, 2021;, btab396, https://doi.org/10.1093/bioinformatics/btab396
