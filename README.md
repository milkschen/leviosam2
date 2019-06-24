# Liftover tool for pairs of genomic sequences

**Purpose**: a C++ library to translate coordinates between two closely related
genomic sequences given a predefined pairwise alignment.

uses succinct data structures as implemented by
[sdsl-lite](https://github.com/simongog/sdsl-lite)

## Usage (Command line)

To serialize liftover information between a reference sequence and given alternative sequence described in a VCF file:

```
$ ./liftover serialize -v <vcf> -s <sample_name> -p <output prefix>
liftover file saved to <output prefix>.lft
```

To lift over coordinates given from a serialized `.lft` file

```
$ ./liftover lift -l <.lft file> -p <output prefix>
lifted coordinates saved to <output prefix>.sam
```

To lift over coordinates on the fly:

```
$ ./liftover lift -v <vcf> -s <sample_name> -p <output prefix>
lifted coordinates saved to <output prefix>.sam
```

## Usage (C++)

To build a liftover between a reference genome and an alternative genotype
given a *VCF* file containing a genotype for a sample over a set of
variants(`FMT/GT` field):

```
#include <cstdio>
#include <htslib/vcf.h>
#include "liftover.hpp"

int main() {
    const char* fname = "data/dna.vcf";
    vcfFile* fp = bcf_open(fname, "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    lift::LiftMap l(fp, hdr, "sample_name");
}
```

To query the equivalent reference position for a given haplotype position:

```
    l.alt_to_ref("contig_name", 8)); // give the contig name and the position on the contig
```

To serialize the liftover to a file:

```
    std::ofstream o("liftover.lft");
    l.serialize(o);
    o.close();
```

To load a liftover from a serialized file

```
    std::ifstream in("liftover.lft");
    lift::Lift l2(in);
    in.close();
```

## Dependencies

- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [htslib](https://github.com/samtools/htslib)

## Currently Supported:

- create a liftover from VCF file w/ FMT/GT field for a specified sample
- reading/writing to SAM files converted the POS field appropriately

## TODO

- ~~convert CIGAR strings for alignments (requires additional SNP information)~~ (additional testing required)
- recalculate read-pair information
- recalculate MAPQ
