# Alignment with alternative reference genomes

In this tutorial, we will demonstrate an alignment pipeline that utilized an alternative reference genome and `levioSAM` to enchance accuracy.
We use 8 threads whenever the program supports multi-threaded processing.

## Software and datasets involved in this pipeline
Software:
- levioSAM
- mason2 simulator (v2.0.9)
- bcftools (v1.11)
- samtools (v1.11)
- Bowtie 2 (v2.4.2) or Bwa-mem (v0.7.17-r1188)

Datasets:
- GRCh38
- Small variant calls (SNVs and indels) from the 1000 Genomes Project
- A chromosome name map file (map chromosome name "Z" to "chrZ")
- A chromosome length map file (provides the length for each chromosome)

```shell
DIR_DATA="resources"
# GRCh38
wget -P ${DIR_DATA}/ ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
bgzip -d ${DIR_DATA}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
samtools faidx ${DIR_DATA}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna chr21 > ${DIR_DATA}/chr21.fa
# chromosome name map file
wget -P ${DIR_DATA}/ https://raw.githubusercontent.com/langmead-lab/reference_flow/master/resources/GRCh38.chrom_map
# chromosome length map file
wget -P ${DIR_DATA}/ https://raw.githubusercontent.com/langmead-lab/reference_flow/master/resources/GRCh38.length_map
# 1KGP genotypes
wget -P ${DIR_DATA}/1kg_vcf/ http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/release/20190312_biallelic_SNV_and_INDEL/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz
```


## Simulate data

We will use a simulated dataset in this tutorial. The advantage of simulated sequencing data is the advent of a ground truth.

We use the genotypes from the first haplotype of individual NA12878 from the 1000 Genomes Project (1KGP).
We only use chr21 data for simplicity.
We simulate 100000 paired-end reads, using the Illumina model.

```sh
SAMPLE="NA12878"
mkdir ${SAMPLE}
bcftools view -s ${SAMPLE} ${DIR_DATA}/1kg_vcf/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | bcftools annotate --rename-chrs ${DIR_DATA}/GRCh38.chrom_map -O z -o ${SAMPLE}/chr21-${SAMPLE}.vcf.gz; bcftools index ${SAMPLE}/chr21-${SAMPLE}.vcf.gz
# Generate NA12878 chr21 hapA. 
# This will be used as the personalized reference genome later.
bcftools consensus -f ${DIR_DATA}/chr21.fa -H 1 -o ${SAMPLE}/${SAMPLE}-chr21-hapA.fa ${SAMPLE}/chr21-${SAMPLE}.vcf.gz
# Simulate reads using mason2
cd ${SAMPLE}/
mason_simulator -ir ${SAMPLE}-chr21-hapA.fa -o ${SAMPLE}-chr21-hapA.1.fq -or ${SAMPLE}-chr21-hapA.2.fq -oa ${SAMPLE}-chr21-hapA-mason.sam -n 100000 --num-threads 8
cd ../
```



## Build major-allele reference

In a major-allele reference, minor alleles in the original reference genome (here GRCh38) are replaced by major alleles cateloged in a reference panel. The resulting reference can be interpreted as the _centroid_ reference genome -- on average it has a closer distance to individuals in the population.

```sh
bcftools view -O z -q 0.5000001 -G ${DIR_DATA}/1kg_vcf/ALL.chr21.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz | bcftools annotate --rename-chrs ${DIR_DATA}/GRCh38.chrom_map -O z -o chr21-major.vcf.gz 
bcftools index chr21-major.vcf.gz
bcftools consensus -f ${DIR_DATA}/chr21.fa -o chr21-major.fa chr21-major.vcf.gz
```



## Prepare levioSAM indexes

```sh
leviosam serialize -v chr21-major.vcf.gz -k ${DIR_DATA}/GRCh38.length_map -p chr21-major
leviosam serialize -v ${SAMPLE}/chr21-${SAMPLE}.vcf.gz -k ${DIR_DATA}/GRCh38.length_map -p ${SAMPLE}/chr21-${SAMPLE} -s ${SAMPLE} -g 0
```



## Align with Bowtie 2

When using major-allele reference and the personalized reference, we pipe the alignment results directory to `leviosam` for smaller disk space usage.
If you'd like to keep the original alignment. You can use typical alignment commands, and use the `leviosam lift -l <lft> -a <input> -p <out_prefix>` command.

```sh
# Build Bowtie 2 indexes
mkdir bt2_indexes
bowtie2-build --threads 8 ${DIR_DATA}/chr21.fa bt2_indexes/chr21
bowtie2-build --threads 8 chr21-major.fa bt2_indexes/chr21-major
bowtie2-build --threads 8 ${SAMPLE}/${SAMPLE}-chr21-hapA.fa bt2_indexes/${SAMPLE}-chr21-hapA
# Alignment
mkdir bt2_aln
R1="${SAMPLE}/${SAMPLE}-chr21-hapA.1.fq"
R2="${SAMPLE}/${SAMPLE}-chr21-hapA.2.fq"
# GRCh38
bowtie2 -p 8 -x bt2_indexes/chr21 -1 ${R1} -2 ${R2} -S bt2_aln/${SAMPLE}-chr21-grch38-bt2.sam
# Major
bowtie2 -p 8 -x bt2_indexes/chr21-major -1 ${R1} -2 ${R2} | leviosam lift -l chr21-major.lft -t 8 -p bt2_aln/${SAMPLE}-chr21-major-bt2
# Personalized-hapA
bowtie2 -p 8 -x bt2_indexes/${SAMPLE}-chr21-hapA -1 ${R1} -2 ${R2} | leviosam lift -l ${SAMPLE}/chr21-${SAMPLE}.lft -t 8 -p bt2_aln/${SAMPLE}-chr21-per-bt2
```

## Align with bwa-mem

Similar to the Bowtie 2 part, we pipe alignment results directory to `leviosam` if liftover is needed.

```sh
# Build bwa mem indexes
mkdir bwa_indexes
bwa index -p bwa_indexes/chr21.fa ${DIR_DATA}/chr21.fa
bwa index -p bwa_indexes/chr21-major.fa chr21-major.fa
bwa index -p bwa_indexes/${SAMPLE}-chr21-hapA.fa ${SAMPLE}/${SAMPLE}-chr21-hapA.fa
# Alignment
mkdir bwa_aln
R1="${SAMPLE}/${SAMPLE}-chr21-hapA.1.fq"
R2="${SAMPLE}/${SAMPLE}-chr21-hapA.2.fq"
# GRCh38
bwa mem -t 8 bwa_indexes/chr21.fa ${R1} ${R2} > bwa_aln/${SAMPLE}-chr21-grch38-bwa.sam
# Major
bwa mem -t 8 bwa_indexes/chr21-major.fa ${R1} ${R2} | leviosam lift -l chr21-major.lft -t 8 -p bwa_aln/${SAMPLE}-chr21-major-bwa
# Personalized-hapA
bwa mem -t 8 bwa_indexes/${SAMPLE}-chr21-hapA.fa ${R1} ${R2} | leviosam lift -l ${SAMPLE}/chr21-${SAMPLE}.lft -t 8 -p bwa_aln/${SAMPLE}-chr21-per-bwa
```

## Evaluate

To evaluate alignment concordance, we first lift the simulation profile (in SAM format) so that all results use GRCh38 coordinates.

```sh
leviosam lift -a ${SAMPLE}/${SAMPLE}-chr21-hapA-mason.sam -l ${SAMPLE}/chr21-${SAMPLE}.lft -p ${SAMPLE}/${SAMPLE}-chr21-hapA-mason-lifted -t 8
```

We use the `compare_sam.py` [script](https://github.com/alshai/levioSAM/blob/master/scripts/compare_sam.py) to evaluate alignment concordance. First set `$PATH_LVSAM` to the levioSAM directory, and then we're good to start run evaluation:

```sh
# Bowtie 2
for fn in bt2_aln/${SAMPLE}-chr21-grch38-bt2.sam bt2_aln/${SAMPLE}-chr21-major-bt2.sam bt2_aln/${SAMPLE}-chr21-per-bt2.sam; do
python $PATH_LVSAM/scripts/compare_sam.py -b ${SAMPLE}/${SAMPLE}-chr21-hapA-mason-lifted.sam -q ${fn} -o ${fn}.cmp_rpt -p 10
samtools view -f 4 ${fn} | wc -l > ${fn}.unmapped_count
done
# bwa-mem
for fn in bwa_aln/${SAMPLE}-chr21-grch38-bwa.sam bwa_aln/${SAMPLE}-chr21-major-bwa.sam bwa_aln/${SAMPLE}-chr21-per-bwa.sam; do
python $PATH_LVSAM/scripts/compare_sam.py -b ${SAMPLE}/${SAMPLE}-chr21-hapA-mason-lifted.sam -q ${fn} -o ${fn}.cmp_rpt -p 10
samtools view -f 4 ${fn} | wc -l > ${fn}.unmapped_count
done
```

### Bowtie 2 results

```sh
for i in `ls bt2_aln/*.cmp_rpt`; do ls $i; cat $i | grep -A 2 'Position'; done
for i in `ls bt2_aln/*count`; do ls $i; cat $i; done
```

|              | # unaligned reads | % correct reads | % correct alignment |
| ------------ | ----------------- | --------------- | ------------------- |
| GRCh38       | 26                | 0.955015        | 0.955139            |
| Major        | 10                | 0.955195        | 0.955243            |
| Personalized | 0                 | 0.955335        | 0.955335            |

### Bwa mem results

```sh
for i in `ls bwa_aln/*.cmp_rpt`; do ls $i; cat $i | grep -A 2 'Position'; done
for i in `ls bwa_aln/*count`; do ls $i; cat $i; done
```

|              | # unaligned reads | % correct reads | % correct alignment |
| ------------ | ----------------- | --------------- | ------------------- |
| GRCh38       | 0                 | 0.954025        | 0.954025            |
| Major        | 0                 | 0.954105        | 0.954105            |
| Personalized | 0                 | 0.954305        | 0.954305            |

We show that in using major-allele reference improves alignment accuracy compared to using GRCh38. If a personalized reference genome is available, the accuracy can be further boosted.
