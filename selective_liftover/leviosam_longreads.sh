# LevioSAM pipeline to lift and re-align long reads
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
DEFER_DEST_BED=""
COMMIT_SOURCE_BED=""
ALLOWED_GAPS=500
ALN=minimap2


while getopts a:C:D:g:F:i:L:M:o:r:R:t:T: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFER_DEST_BED=" -D ${OPTARG}";;
        R) COMMIT_SOURCE_BED=" -r ${OPTARG}";;
        F) REF=${OPTARG};;
        g) ALLOWED_GAPS=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
        o) PFX=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Aligner: ${ALN}";
echo "Reference: ${REF}";
echo "Aligner read group: ${ALN_RG}";
echo "LevioSAM software: ${LEVIOSAM}";
echo "LevioSAM index: ${CLFT}";
echo "Allowed gaps: ${ALLOWED_GAPS}";
echo "BED where reads get deferred: ${DEFER_DEST_BED}";
echo "BED where reads get discarded: ${COMMIT_SOURCE_BED}";
echo "Num. threads: ${THR}";

if [[ ! ${ALN} =~ ^(minimap2|winnowmap2)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: minimap2, winnowmap2"
    exit
fi

# Lifting over using leviosam
if [ ! -s ${PFX}-committed.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log \
            ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
            -S lifted -G ${ALLOWED_GAPS} \
            ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
    else
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
        -S lifted -G ${ALLOWED_GAPS} \
        ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
    fi
fi

# Convert deferred reads to FASTQ
if [ ! -s ${PFX}-deferred.fq.gz ]; then
    samtools fastq ${PFX}-deferred.bam | bgzip > ${PFX}-deferred.fq.gz
fi

# Re-align deferred reads
if [ ! -s ${PFX}-realigned.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o aln_deferred.time_log ${ALN} -ax asm20 -t ${THR} ${REF} ${PFX}-deferred.fq.gz | samtools view -hbo ${PFX}-realigned.bam
    else
        ${ALN} -ax asm20 -t ${THR} ${REF} ${PFX}-deferred.fq.gz | samtools view -hbo ${PFX}-realigned.bam
    fi
fi

# Merge and sort
if [ ! -s ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge_and_sort.time_log \
            samtools cat ${PFX}-committed.bam ${PFX}-realigned.bam | samtools sort -@ ${THR} -o ${PFX}-final.bam
    else
        samtools cat ${PFX}-committed.bam ${PFX}-realigned.bam | samtools sort -@ ${THR} -o ${PFX}-final.bam
    fi
fi

