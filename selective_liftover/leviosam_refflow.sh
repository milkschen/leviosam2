set -x

INPUT="HG002.bt2.chm13_v1.0.bam"
PFX="HG002.bt2.chm13_v1.0_to_h38"
TH=8
STH=4
BT2_IDX="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/fasta/grch38/bt2/GCA_000001405.15_GRCh38_no_alt_analysis_set"
BT2_RG="--rg-id chm13to38 --rg SM:HG002 --rg LB:chm13 --rg PL:ILLUMINA --rg DS:novaseq --rg PU:novaseq"
LEVIOSAM="$NC2/levioSAM/leviosam"
CLFT="$NC2/fasta/t2t-chm13-v1.0.hg38.clft"
TIME="/home-1/cnaechy1@jhu.edu/bin/time-1.9"

# Check inputs
if [ ! -f ${PFX}.bam ]; then
    ${TIME} -v -o lift.time_log ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${TH} -p ${PFX} -O bam -S mapq
fi

if [ ! -f ${PFX}-deferred-collate.bam ]; then
    ${TIME} -v -o collate.time_log samtools collate -@ ${TH} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
fi
FQSIZE=$(stat -c%s "${PFX}-deferred-R1.fq")
if (( ${FQSIZE} > 0 )); then
    echo "Skip converting FASTQ"
else
    ${TIME} -v -o to_fastq.time_log samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
fi

if [ ! -f ${PFX}-deferred-re_aligned-paired.bam ]; then
    ${TIME} -v -o aln_paired.time_log bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
fi
if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    ${TIME} -v -o aln_unpaired.time_log bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
fi

if [ ! -f ${PFX}-merged.bam ]; then
    ${TIME} -v -o merge.time_log samtools merge ${PFX}-merged.bam ${PFX}.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam
fi

if [ ! -f ${PFX}-final.bam ]; then
    ${TIME} -v -o sort_all.time_log samtools sort -@ ${TH} -o ${PFX}-final.bam ${PFX}-merged.bam
fi
