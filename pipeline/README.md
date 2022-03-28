# Selective liftover pipelines (ChainMap)

Combining leviosam with a selective re-alignment strategy can improve accuracy.

## Short read (paired-end) pipeline (`leviosam_shortreads_pe.sh`)

This pipeline is designed for paired-end short reads (Illumina). 
We first lift all reads using `leviosam lift`. 
Reads are split into _committed_ and _deferred_ groups using a set of filtration criteria. 
We then use `leviosam collate` to make sure all reads are properly paired. 
If one pair is assigned as both committed and deferred, we'll put both segments into the deferred group.
The deferred reads are then converted to FASTQ and re-aligned.
Committed and re-aligned deferred reads are concatenated and sorted as the final output.

* Supported aligners: [bwa mem](https://github.com/lh3/bwa) and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 

Example:
```
sh /home/cnaechy1/data_blangme2/naechyun/levioSAM/selective_liftover/leviosam_refflow_mappability.sh \
    -i pe_short_reads.bam \
    -o lifted_pe_short_reads \
    -b GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -C hg19ToHg38.over.clft \
    -t 16 \
    -D hg19-grch38-lt98-map_reduction_k100_e0.01-merged.bed
```


## Long read pipeline (`leviosam_longreads.sh`)

The long read pipeline is similar to the short read pipeline, with the following differences:

- Reads are single-end
- Supported aligners: [minimap2](https://github.com/lh3/minimap2) and [winnowmap2](https://github.com/marbl/Winnowmap)

Example:
```
sh leviosam_longreads.sh \
    -i long_reads.bam \
    -o lifted_long_reads \
    -F GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -C hg19ToHg38.over.clft \
    -t 16
```
