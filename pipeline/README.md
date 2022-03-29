# LevioSAM2 pipelines

## Dependencies

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
    -a minimap2 -g 1000 -H 100 -x ../configs/pacbio_all.yaml \
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

# Download FASTAs and pre-built indexes
## GRCh38
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
## Pre-build GRCh38 indexes
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
## CHM13v2
wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fasta
bowtie2-build --threads 16 chm13v2.0/chm13v2.0 chm13v2.0.fasta

# Download levioSAM2 resources

# Align to CHM13
bowtie2 -p 16 -x GRCh38_noalt_as/GRCh38_noalt_as -1 {} -2 {} | samtools view -hbo {}.bam

# Run the levioSAM2 pipeline
```

