# LevioSAM2 workflow

> [!WARNING]
> **CLI Flag Differences between Workflows and the C++ Binary:**
> The flag names in these workflow scripts (`leviosam2.py` / `leviosam2.sh`) differ from those of the underlying C++ binary (`leviosam2 lift`). If you run the C++ binary directly, verify the option mappings:
> * **Gap size threshold:** Use `-g` in Python/Bash scripts, but **`-G`** in the C++ binary. (The C++ binary's `-g` is reserved for VCF haplotype).
> * **Sequence / Aligner / Sample flags:** `-s` specifies sequence type in Python, source index path in Bash, and VCF sample name in the C++ binary.
> * **Commit / Suppress annotations:** Use `--lift_bed_commit_source` (Python) or `-R` (Bash), but **`-r`** in the C++ binary.
> * **Defer rules (`-m` and `-p`):** In Bash, `-m` is the template size cutoff and `-p` is the clipped fraction cutoff. In the C++ binary, `-m` enables MD/NM tags, and `-p` is the **output prefix**.
> Check `leviosam2 lift -h` for details.

## Dependencies

All dependent software is included in our Docker/Singularity container.
Users of other installation approaches may install the dependencies using
conda, pre-built binary, or building from source.

- Aligner: [minimap2](https://github.com/lh3/minimap2)/[winnowmap2](https://github.com/marbl/Winnowmap)/[bwa mem](https://github.com/lh3/bwa)/[Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- samtools, bgzip

## Paired-end workflows

- Supported aligners: [bwa mem](https://github.com/lh3/bwa) and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

The annotation options in the examples are optional: `--lift_bed_commit_source` / `--lift_bed_defer_target` in Python, or `-R` / `-D` in Bash. Omit them if you do not have annotation BED files.

Bowtie 2 example:

```shell
python leviosam2.py \
    -t 16 \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -C chm13v2-grch38.clft \
    -f grch38.fna \
    -fi bt2/grch38.fna \
    --sequence_type ilmn_pe \
    -a bowtie2 \
    --use_preset \
    --lift_bed_commit_source suppress_annotations.bed \
    --lift_bed_defer_target defer_annotations.bed
```

Alternative with the Bash workflow (not recommended):

```shell
bash leviosam2.sh \
    -a bowtie2 -A -10 -q 10 -H 5 \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -f grch38.fna \
    -b bt2/grch38 \
    -C chm13v2-grch38.clft \
    -t 16 \
    -D defer_annotations.bed \
    -R suppress_annotations.bed
```

BWA MEM example:

```shell
python leviosam2.py \
    -t 16 \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -C chm13v2-grch38.clft \
    -f grch38.fna \
    -fi bwa/grch38.fna \
    --sequence_type ilmn_pe \
    -a bwamem \
    --use_preset \
    --lift_bed_commit_source suppress_annotations.bed \
    --lift_bed_defer_target defer_annotations.bed
```

Alternative with the Bash workflow (not recommended):

```shell
bash leviosam2.sh \
    -a bwamem -A 100 -q 30 -H 5 \
    -i ilmn-pe.bam \
    -o ilmn-pe-lifted \
    -f grch38.fna \
    -b bwa/grch38.fna \
    -C chm13v2-grch38.clft \
    -t 16 \
    -D defer_annotations.bed \
    -R suppress_annotations.bed
```

## Single-end workflows

For both short and long reads. Different parameters are recommended for each sequence type.

- Supported aligners: [minimap2](https://github.com/lh3/minimap2), [winnowmap2](https://github.com/marbl/Winnowmap),
[bwa mem](https://github.com/lh3/bwa) and [Bowtie 2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)

- `-g` sets the max allowed size of overlapping chain gaps of an alignment, `-H` is the max edit distance of an alignment (post-realignment) to be committed. Please adjust these parameters according to your long read platform and library preparation

Minimap2 example with Pacbio:

```shell
python leviosam2.py \
    -t 16 \
    -i pacbio.bam \
    -o pacbio-lifted \
    -C chm13v2-grch38.clft \
    -f grch38.fna \
    --sequence_type pb_hifi \
    -a minimap2 \
    --use_preset \
    --lift_realign_config ../configs/pacbio_all.yaml \
    --lift_bed_commit_source suppress_annotations.bed
```

Alternative with the Bash workflow (not recommended):

```shell
bash leviosam2.sh \
    -a minimap2 -g 1000 -H 100 -S -x ../configs/pacbio_all.yaml \
    -l map-hifi \
    -i pacbio.bam \
    -o pacbio-lifted \
    -f grch38.fna \
    -C chm13v2-grch38.clft \
    -t 16 \
    -R suppress_annotations.bed
```

Minimap2 example with Nanopore (ultra-long reads are not fully supported):

```shell
python leviosam2.py \
    -t 16 \
    -i ont.bam \
    -o ont-lifted \
    -C chm13v2-grch38.clft \
    -f grch38.fna \
    --sequence_type ont \
    -a minimap2 \
    --use_preset \
    --lift_realign_config ../configs/ont_all.yaml \
    --lift_bed_commit_source suppress_annotations.bed
```

Alternative with the Bash workflow (not recommended):

```shell
bash leviosam2.sh \
    -a minimap2 -g 1500 -H 6000 -S -x ../configs/ont_all.yaml \
    -l map-ont \
    -i ont.bam \
    -o ont-lifted \
    -f grch38.fna \
    -C chm13v2-grch38.clft \
    -t 16 \
    -R suppress_annotations.bed
```

## Output files

The Python workflow writes an indexed BAM by default. Use `-O cram` to write an indexed CRAM instead.

| File | Contents |
| --- | --- |
| `<out_prefix>-final.bam` and `.bam.bai` | Final lifted and realigned reads in target coordinates (default). |
| `<out_prefix>-final.cram` and `.cram.crai` | Final reads in target coordinates when using `-O cram`. |
| `<out_prefix>-unliftable.bam` | Deferred and unliftable reads in source coordinates. See [#20](https://github.com/milkschen/leviosam2/issues/20). |

The workflow accepts BAM and CRAM output. For SAM output in quick experiments, run `leviosam2 lift -O sam` directly.

Temporary files are removed by default. Use `--keep_tmp_files` in Python or `-K` in Bash to retain intermediate alignments and deferred FASTQ files. The Bash workflow produces BAM output only.

## Case Study

Please make sure `leviosam2` and `leviosam2.sh` is in your `$PATH` with command `leviosam2 -h` and `sh leviosam2.sh -h`.
If it fails, please use the following command to set `$PATH` (paths need to be adjusted locally) and make sure you see the help message.

```shell
PATH=<path/to/leviosam2/build>:<path/to/leviosam2/workflow>:${PATH}
```

```shell
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
bowtie2 -p 16 -x chm13v2.0/chm13v2.0 -1 HG002.novaseq.pcr-free.0_3x-R1.fq.gz -2 HG002.novaseq.pcr-free.0_3x-R2.fq.gz | samtools sort -@ 4 -o ilmn-pe-chm13v2.bam

# Run the levioSAM2 workflow
leviosam2 index -c chm13v2-grch38/chm13v2-grch38.chain -p chm13v2-grch38/chm13v2-grch38 -F GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.fai
bash ../leviosam2.sh \
    -a bowtie2 -A -10 -q 10 -H 5 \
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

### Read groups

For Bowtie2, pass its read-group options together in `--read_group`, for example
`--read_group "--rg-id foo --rg 'SM:sample name'"`. Quote any option value
containing spaces within that string. For BWA, minimap2, winnowmap2, and
strobealign, pass the read-group value as a single string.
