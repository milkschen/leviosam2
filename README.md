# levioSAM lifts variant-aware alignments to the reference genome

Use a VCF file containing alternative haplotype information to lift SAM/BAM alignments
to the reference sequence


## Supported Features:

- Serialization of VCF file for faster parsing

- Converting SAM/BAM records from haplotype to reference, including:
    - CIGAR string
    - MD:Z and NM:i tag
    - paired read information

- Multithreading support

## Building and installation

The easiest way to install levioSAM is by using [Conda](https://docs.conda.io/en/latest/):

```
conda install -c conda-forge -c bioconda leviosam
```

We support a variety of other methods for getting levioSAM to work, including Docker, CMake and Make. See
[INSTALL.md](INSTALL.md) for more details.

## Usage (command line)

We highly suggest you normalize and left-align your VCF file before using it to lift your alignments:
```
bcftools norm <VCF> > <output VCF>
```

Run this command to lift your alignments:
```
$ ./levioSAM lift -t <nthreads> -a <sam> -v <vcf> -s <sample_name> -g <haplotype (0 or 1)> -p <output prefix>
```

`<sample_name>` is the name of the individual in your VCF whose genotype you want to use for the lift-over.
if `-s` is not specified, then it will be assumed that all the alternate alleles in all VCF are present in the variant-aware
reference.

The lifted coordinates will be saved in SAM format to `<output prefix>.sam`. If `-p` is not specified then the SAM file
will be output to `stdout`.

You can speed up the lift-over by serializing the lift-over structure beforehand to avoid re-parsing the entire VCF. This
is useful if you want to use the same VCF to run the command multiple times.
```
$ ./levioSAM serialize -v <vcf> -s <sample_name> -p <output prefix>
```
The levioSAM file will saved to `<output prefix>.lft`.

You can then pass the serialized structure into lift-over by using the `-l` option instead of the `-v` option:
```
$ ./levioSAM lift -a <sam> -l <lft> -p <output prefix>
```

A list of common options:

- To read from stdin, use `-a -` or exclude `-a`.
- To use multiple threads, use `-t <threads>`.
- To lift `NM:i` and `MD:z` tags, add `-m` (this will be slow if your alignments are not sorted).
- To write output as a BAM file, use `-O bam` (the lifted file will be `<output prefix>.bam`).

## Example (Command line)

The `testdata` directory contains some toy data that you can play around with. It also
contains a `.lft` file that will let you lift alignments from the [major-allele reference](https://doi.org/10.1371/journal.pgen.1002280) to GRCh38.
Though we recommend you build your own `.lft` file from a VCF for analyses.

Here's an example using a small set of Bowtie 2 paired-end alignments against the major-allele reference

```
leviosam -a testdata/bt2-paired_end-major.sam -l testdata/wg-maj.lft > out.sam
```

We provide instructions of how to use levioSAM in common variant-aware reference pipelines (major-allele reference and personalized reference) in the [levioSAM wiki](https://github.com/alshai/levioSAM/wiki/Alignment-with-variant-aware-reference-genomes).

## Example (C++)

Use the `LiftMap` class to generate lift-over information between a reference genome and an alternative genotype.

```
#include <cstdio>
#include <tuple>
#include <unordered_map>
#include <htslib/vcf.h>
#include <vector>
#include "levioSAM.hpp"

int main() {
    const char* fname = "data/dna.vcf";
    vcfFile* fp = bcf_open(fname, "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    std::string sample_name =  ""; // GT=1 will be assumed for all variants if left blank
    std::string haplotype = "0"; // set this to "0" or "1" if sample_name is not blank
    // populate name_map if the contig names in the VCF file do not match that of the genome
    // (e.g.. chr1 vs 1).
    std::vector<std::pair<std::string,std::string>> name_map;
    // populate contig_lengths if the chromsome lengths are not specified in the VCF file.
    std::unordered_map<std::string> contig_lengths;
    lift::LiftMap l(fp, hdr, sample_name, haplotype, name_map, contig_length);
}
```

To query the equivalent reference position for a given haplotype position:

```
    std::vector<std::string> c; // for compatibility with multi-thread
    std::mutex m; // for compatibility with multi-thread
    l.s2_to_s1("contig_name", 8, &c, &m)); // give the contig name and the position on the contig
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

## Publications

Mun T., Chen N., Langmead B. ["LevioSAM: Fast lift-over of alternate reference alignments"](https://doi.org/10.1101/2021.02.05.429867). [BioRxiv](https://www.biorxiv.org/). 2021
