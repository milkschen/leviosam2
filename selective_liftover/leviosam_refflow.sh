set -x

INPUT="HG002.bt2.chm13_v1.0.bam"
PFX="HG002.bt2.chm13_v1.0_to_h38"
BT2_IDX="/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/fasta/grch38/bt2/GCA_000001405.15_GRCh38_no_alt_analysis_set"
BT2_RG="--rg-id chm13to38 --rg SM:HG002 --rg LB:chm13 --rg PL:ILLUMINA --rg DS:novaseq --rg PU:novaseq"
CLFT="$NC2/fasta/t2t-chm13-v1.0.hg38.clft"
THR=8
LEVIOSAM="$NC2/levioSAM/leviosam" # leviosam
TIME="/home-1/cnaechy1@jhu.edu/bin/time-1.9" # time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step

while getopts i:o:b:r:C:t flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        o) PFX=${OPTARG};;
        b) BT2_IDX=${OPTARG};;
        r) BT2_RG=${OPTARG};;
        C) CLFT=${OPTARG};;
        t) THR=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Bowtie2 indexes prefix: ${BT2_IDX}";
echo "Bowtie2 read group: ${BT2_RG}";
echo "LevioSAM index: ${CLFT}";
TH=$(( ${THR} * 2/3 ))
STH=$(( ${THR} - ${TH} ))
echo "Num. threads: ${THR} (${TH}/${STH})";

# Check inputs
if [ ! -f ${PFX}.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq
    else
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq
    fi
fi

if [ ! -f ${PFX}-deferred-collate.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o collate.time_log samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    else
        samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    fi
fi
FQSIZE=$(stat -c%s "${PFX}-deferred-R1.fq")
if (( ${FQSIZE} > 0 )); then
    echo "Skip converting FASTQ"
else
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o to_fastq.time_log samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    else
        samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    fi
fi

if [ ! -f ${PFX}-deferred-re_aligned-paired.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o aln_paired.time_log bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
    else
        bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
    fi
fi
if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o aln_unpaired.time_log bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
    else
        bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
    fi
fi

if [ ! -f ${PFX}-merged.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge.time_log samtools merge ${PFX}-merged.bam ${PFX}.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam
    else
        samtools merge ${PFX}-merged.bam ${PFX}.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam
    fi
fi

if [ ! -f ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o sort_all.time_log samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    else
        samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    fi
fi
