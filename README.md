# Liftover tool for pairs of genomic sequences

**Purpose**: a C++ library to translate coordinates between two closely related
genomic sequences given a predefined pairwise alignment.

uses succinct data structures as implemented by
[sdsl-lite](https://github.com/simongog/sdsl-lite)


## Usage

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

# Dependencies

- [sdsl-lite](https://github.com/simongog/sdsl-lite)
- [htslib](https://github.com/samtools/htslib)
