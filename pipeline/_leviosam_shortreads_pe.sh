# LevioSAM2 pipeline to lift and re-align paired-end short reads
#
# Authors: Nae-Chyun Chen
#
# Distributed under the MIT license
# https://github.com/milkschen/leviosam2
#

ALN_RG=""
THR=$(nproc)
LEVIOSAM=leviosam2
TIME=time # GNU time
MEASURE_TIME=0 # Set to a >0 value to measure time for each step
ALLOWED_GAPS=0
FRAC_CLIPPED=0.95
ISIZE=1000
HDIST=5
DEFER_DEST_BED=""
COMMIT_SOURCE_BED=""
REALN_CONFIG=""

# Default parameters for Bowtie 2 (local)
# ALN=bowtie2
# MAPQ=10
# ALN_SCORE=200

# Default parameters for BWA-MEM
ALN=bwamem
MAPQ=30
ALN_SCORE=100


while getopts a:A:b:C:D:f:H:g:i:L:M:o:q:r:R:t:T:x: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        A) ALN_SCORE=${OPTARG};;
        b) ALN_IDX=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFER_DEST_BED=" -D ${OPTARG}";;
        f) REF=${OPTARG};;
        H) HDIST=${OPTARG};;
        g) ALLOWED_GAPS=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
        o) PFX=${OPTARG};;
        q) MAPQ=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        R) COMMIT_SOURCE_BED=" -r ${OPTARG}";;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
        x) REALN_CONFIG=" -x ${OPTARG}";;
    esac
done

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

TT=""
if (( ${MEASURE_TIME} > 0 )); then
    TT="${TIME} -v -ao leviosam.time_log "
fi

if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit 1
fi

set -xp

# Lifting over using leviosam
if [ ! -s ${PFX}-committed.bam ]; then
    ${TT} ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
    -S mapq:${MAPQ} -S isize:${ISIZE} -S aln_score:${ALN_SCORE} \
    -S clipped_frac:${FRAC_CLIPPED} -S hdist:${HDIST} -G ${ALLOWED_GAPS} \
    ${REALN_CONFIG} \
    -m -f ${REF} ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
fi

# Collate
if [ ! -s ${PFX}-paired-deferred-R1.fq.gz ]; then
    ${TT} ${LEVIOSAM} collate -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
fi

# Re-align deferred reads
if [ ! -s ${PFX}-paired-realigned.bam ]; then
    if [[ ${ALN} == "bowtie2" ]]; then
        ${TT} bowtie2 --local ${ALN_RG} -p ${THR} -x ${ALN_IDX} \
        -1 ${PFX}-paired-deferred-R1.fq.gz -2 ${PFX}-paired-deferred-R2.fq.gz | \
        ${TT} samtools view -hb > ${PFX}-paired-realigned.bam
    else
        ${TT} bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} \
        ${PFX}-paired-deferred-R1.fq.gz ${PFX}-paired-deferred-R2.fq.gz | \
        ${TT} samtools view -hb > ${PFX}-paired-realigned.bam
    fi
fi

# Reference flow-style merging
if [ ! -s ${PFX}-paired-realigned-sorted_n.bam ]; then
    ${TT} samtools sort -@ ${THR} -n -o ${PFX}-paired-realigned-sorted_n.bam ${PFX}-paired-realigned.bam
fi
if [ ! -s ${PFX}-paired-deferred-sorted_n.bam ]; then
    ${TT} samtools sort -@ ${THR} -n -o ${PFX}-paired-deferred-sorted_n.bam ${PFX}-paired-deferred.bam
fi
if [ ! -s ${PFX}-paired-deferred-reconciled.bam ]; then
    ${TT} ${LEVIOSAM} reconcile -s source:${PFX}-paired-deferred-sorted_n.bam -s target:${PFX}-paired-realigned-sorted_n.bam \
        -m -o ${PFX}-paired-deferred-reconciled.bam
fi

# Merge, sort, and clean
if [ ! -s ${PFX}-final.bam ]; then
    ${TT} samtools cat ${PFX}-paired-committed.bam ${PFX}-paired-deferred-reconciled.bam | \
        ${TT} samtools sort -@ ${THR} -o ${PFX}-final.bam
    rm ${PFX}-paired-deferred.bam ${PFX}-paired-deferred-sorted_n.bam 
    rm ${PFX}-paired-realigned.bam ${PFX}-paired-realigned-sorted_n.bam
    rm ${PFX}-paired-deferred-reconciled.bam
    rm ${PFX}-paired-committed.bam ${PFX}-deferred.bam
fi
