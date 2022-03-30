# LevioSAM2 pipelines

## Dependencies

All dependent software is included in our Docker/Singularity container.
Users of other installation approaches may install the dependencies using conda, pre-built binary, or building from source.

- Aligner: [minimap2](https://github.com/lh3/minimap2)/[winnowmap2](https://github.com/marbl/Winnowmap)/[bwa mem](https://github.com/lh3/bwa)/[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 
- samtools, bgzip

## Paired-end pipelines

* Supported aligners: [bwa mem](https://github.com/lh3/bwa) and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 

Bowtie 2 example:
```
sh leviosam2.sh \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -f grch38.fna \
    -b bt2/grch38 \
    -C chm13v2-grch38.clft \
    -t 16
```

BWA MEM example:
```
sh leviosam2.sh \
    -a bwamem -A 100 -q 30 \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -f grch38.fna \
    -b bwa/grch38.fna \
    -C chm13v2-grch38.clft \
    -t 16
```

## Single-end pipelines

For both short and long reads. Different parameters are recommended for each sequence type.

- Supported aligners: [minimap2](https://github.com/lh3/minimap2), [winnowmap2](https://github.com/marbl/Winnowmap),
[bwa mem](https://github.com/lh3/bwa) and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml) 

Minimap2 example:
```
sh leviosam2.sh \
    -a minimap2 -g 1000 -H 100 -S -x ../configs/pacbio_all.yaml \
    -i pacbio.bam \
    -o pacbio-lifted \
    -f grch38.fna \
    -C chm13v2-grch38.clft \
    -t 16
```


## Case Study
```
mkdir case_study
cd case_study

# Download FASTQs
wget https://genome-idx.s3.amazonaws.com/lev/HG002.novaseq.pcr-free.0_3x-R1.fq.gz
wget https://genome-idx.s3.amazonaws.com/lev/HG002.novaseq.pcr-free.0_3x-R2.fq.gz

# Download FASTAs and pre-built indexes
## GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
## Pre-build GRCh38 indexes
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip

## CHM13v2
wget https://genome-idx.s3.amazonaws.com/bt/chm13v2.0.zip
unzip chm13v2.0.zip

# Download levioSAM2 resources
wget https://genome-idx.s3.amazonaws.com/lev/chm13v2-grch38.tar.gz
tar xfzv chm13v2-grch38.tar.gz

# Align to CHM13
bowtie2 -p 16 -x chm13v2.0/chm13v2.0 -1 HG002.novaseq.pcr-free.0_3x-R1.fq.gz -2 HG002.novaseq.pcr-free.0_3x-R2.fq.gz | samtools view -hbo ilmn-pe-chm13v2.bam

# Run the levioSAM2 pipeline
leviosam2 index -c chm13v2-grch38/chm13v2-grch38.chain -p chm13v2-grch38/chm13v2-grch38 -F GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
sh leviosam2.sh \
    -i ilmn-pe-chm13v2.bam \
    -o ilmn-pe-chm13v2_grch38 \
    -f GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
    -b GRCh38_noalt_as/GRCh38_noalt_as \
    -C chm13v2-grch38/chm13v2-grch38.clft \
    -D chm13v2-grch38/chm13v2-grch38-lt0.98-map_reduction_k100_e0.01.bed \
    -R chm13v2-grch38/chm13v2-grch38-source-unliftable-s_5000-winnowmap2-exc_1.bed \
    -t 16
```

The output is `ilmn-pe-chm13v2_grch38-final.bam`.
