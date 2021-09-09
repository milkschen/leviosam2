set -x

PFX="HG002.novaseq.pcr-free.35x.bt2.chm13_v1.0_to_h38"
TH=8
STH=4
BT2_IDX="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/fasta/grch38/bt2/GCA_000001405.15_GRCh38_no_alt_analysis_set"
BT2_RG="--rg-id chm13to38 --rg SM:HG002 --rg LB:chm13 --rg PL:ILLUMINA --rg DS:novaseq --rg PU:novaseq"

# $NC2/levioSAM/leviosam lift -C $NC2/fasta/t2t-chm13-v1.0.hg38.clft -a HG002.novaseq.pcr-free.35x.bt2.chm13_v1.0.bam -t 4 -p HG002.novaseq.pcr-free.35x.bt2.chm13_v1.0_to_h38 -O bam -S mapq[

# Check inputs
if [ ! -f ${PFX}-deferred.bam ]; then
    echo "Missing input: ${PFX}-deferred.bam"
    exit
fi
if [ ! -f ${PFX}.bam ]; then
    echo "Missing input: ${PFX}.bam"
    exit
fi

# if [ ! -f ${PFX}-deferred-sort_n.bam ]; then
    # samtools sort -n -@ ${TH} -o ${PFX}-deferred-sort_n.bam ${PFX}-deferred.bam
if [ ! -f ${PFX}-deferred-collate.bam ]; then
    samtools collate -@ ${TH} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
fi
FQSIZE=$(stat -c%s "${PFX}-deferred-R1.fq")
# if [ ! -f ${PFX}-deferred-R1.fq ]; then
if (( ${FQSIZE} > 0 )); then
    echo "Skip converting FASTQ"
else
    samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
fi

if [ ! -f ${PFX}-deferred-re_aligned-paired.bam ]; then
    bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
fi
if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
fi


if [ ! -f  ${PFX}-committed.bam ]; then
    samtools addreplacerg -r ID:chm13to38 -r LB:chm13 -r PL:ILLUMINA -r SM:HG002 -r PU:novaseq ${PFX}.bam | samtools sort -@ ${STH} -o ${PFX}-committed.bam
    # DS:novaseq
    # samtools addreplacerg -r "@RG	ID:HG002	LB:grch38	PL:illumina	SM:grch38	PU:novaseq" ${PFX}.bam | samtools sort -@ ${STH} -o ${PFX}-committed.bam
fi

if [ ! -f ${PFX}-final.bam ]; then
    samtools merge ${PFX}-final.bam ${PFX}-committed.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam
fi
