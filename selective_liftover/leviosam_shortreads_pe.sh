# LevioSAM pipeline to lift and re-align paired-end short reads
#
# Authors: Nae-Chyun Chen
#
# Distributed under the MIT license
# https://github.com/alshai/levioSAM
#
set -xp

ALN_RG=""
THR=$(nproc)
LEVIOSAM=leviosam
TIME=time # GNU time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step
FRAC_CLIPPED=0.95
ISIZE=1000
DEFER_DEST_BED=""
COMMIT_SOURCE_BED=""

# Default parameters for Bowtie 2 (local)
# ALN=bowtie2
# MAPQ=10
# ALN_SCORE=200

# Default parameters for BWA-MEM
ALN=bwamem
MAPQ=30
ALN_SCORE=100


while getopts a:A:b:C:D:f:i:L:M:o:q:r:R:t:T: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        A) ALN_SCORE=${OPTARG};;
        b) ALN_IDX=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFER_DEST_BED=" -D ${OPTARG}";;
        f) REF=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
        o) PFX=${OPTARG};;
        q) MAPQ=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        R) COMMIT_SOURCE_BED=" -r ${OPTARG}";;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Aligner: ${ALN}";
echo "Aligner indexes prefix: ${ALN_IDX}";
echo "Aligner read group: ${ALN_RG}";
echo "Targer reference: ${REF}";
echo "LevioSAM software: ${LEVIOSAM}";
echo "LevioSAM index: ${CLFT}";
echo "LevioSAM min MAPQ: ${MAPQ}";
echo "LevioSAM min AS: ${ALN_SCORE}";
echo "BED where reads get deferred: ${DEFER_DEST_BED}";
echo "BED where reads get discarded: ${COMMIT_SOURCE_BED}";
echo "Num. threads: ${THR}";

if [[ ${INPUT} == "" ]]; then
    echo "Input is not set properly"
    exit 1
fi
if [[ ${PFX} == "" ]]; then
    echo "Prefix is not set properly"
    exit 1
fi
if [[ ${REF} == "" ]]; then
    echo "Reference (source) is not set properly"
    exit 1
fi

if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit 1
fi

# Lifting over using leviosam
if [ ! -s ${PFX}-committed.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log \
            ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
            -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0 \
            -m -f ${REF} \
            ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
    else
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
        -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0 \
        -m -f ${REF} \
        ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
    fi
fi

# Collate
if [ ! -s ${PFX}-paired-deferred-R1.fq.gz ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o collate.time_log \
            ${LEVIOSAM} collate -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
    else
        ${LEVIOSAM} collate -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
    fi
fi

# Re-align deferred reads
if [ ! -s ${PFX}-paired-realigned.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${TIME} -v -o aln_paired.time_log \
                bowtie2 --local ${ALN_RG} -p ${THR} -x ${ALN_IDX} \
                -1 ${PFX}-paired-deferred-R1.fq.gz -2 ${PFX}-paired-deferred-R2.fq.gz | \
                samtools view -hb > ${PFX}-paired-realigned.bam
        else
            ${TIME} -v -o aln_paired.time_log \
                bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} \
                ${PFX}-paired-deferred-R1.fq.gz ${PFX}-paired-deferred-R2.fq.gz | \
                samtools view -hb > ${PFX}-paired-realigned.bam
        fi
    else
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 --local ${ALN_RG} -p ${THR} -x ${ALN_IDX} \
            -1 ${PFX}-paired-deferred-R1.fq.gz -2 ${PFX}-paired-deferred-R2.fq.gz | \
            samtools view -hb > ${PFX}-paired-realigned.bam
        else
            bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} \
            ${PFX}-paired-deferred-R1.fq.gz ${PFX}-paired-deferred-R2.fq.gz | \
            samtools view -hb > ${PFX}-paired-realigned.bam
        fi
    fi
fi

# Reference flow-style merging
if [ ! -s ${PFX}-paired-realigned-sorted_n.bam ]; then
    samtools sort -@ ${THR} -n -o ${PFX}-paired-realigned-sorted_n.bam ${PFX}-paired-realigned.bam
fi
if [ ! -s ${PFX}-paired-deferred-sorted_n.bam ]; then
    samtools sort -@ ${THR} -n -o ${PFX}-paired-deferred-sorted_n.bam ${PFX}-paired-deferred.bam
fi
if [ ! -s ${PFX}-paired-deferred-cherry_picked.bam ]; then
    ${LEVIOSAM} cherry_pick -s source:${PFX}-paired-deferred-sorted_n.bam -s target:${PFX}-paired-realigned-sorted_n.bam \
        -m -o ${PFX}-paired-deferred-cherry_picked.bam
fi

# Merge and sort
if [ ! -s ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge_and_sort.time_log \
            samtools cat ${PFX}-paired-committed.bam ${PFX}-paired-deferred-cherry_picked.bam | samtools sort -@ ${THR} -o ${PFX}-final.bam
    else
        samtools cat ${PFX}-paired-committed.bam ${PFX}-paired-deferred-cherry_picked.bam | samtools sort -@ ${THR} -o ${PFX}-final.bam
    fi
fi
