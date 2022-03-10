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
ALLOWED_GAPS=1000
HDIST="-S hdist:100"
ALN=minimap2
REALN_CONFIG=""


while getopts a:C:D:g:f:H:i:L:M:o:r:R:t:T:x: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFER_DEST_BED=" -D ${OPTARG}";;
        f) REF=${OPTARG};;
        H) HDIST=" -S hdist:${OPTARG}";;
        g) ALLOWED_GAPS=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
        o) PFX=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        R) COMMIT_SOURCE_BED=" -r ${OPTARG}";;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
        x) REALN_CONFIG=" -x ${OPTARG}";;
    esac
done

if [[ ! ${ALN} =~ ^(minimap2|winnowmap2)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: minimap2, winnowmap2"
    exit
fi

TT=""
if (( ${MEASURE_TIME} > 0 )); then
    TT="${TIME} -v -ao leviosam.time_log "
fi

# Lifting over using leviosam
if [ ! -s ${PFX}-committed.bam ]; then
    ${TT} ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
    -S lifted ${HDIST} -G ${ALLOWED_GAPS} \
    ${REALN_CONFIG} \
    -m -f ${REF} \
    ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
fi

# Convert deferred reads to FASTQ
if [ ! -s ${PFX}-deferred.fq.gz ]; then
    ${TT} samtools fastq ${PFX}-deferred.bam | \
    ${TT} bgzip > ${PFX}-deferred.fq.gz
fi

# Re-align deferred reads
if [ ! -s ${PFX}-realigned.bam ]; then
    ${TT} ${ALN} -ax map-hifi --MD -t ${THR} -R ${ALN_RG} ${REF} ${PFX}-deferred.fq.gz | \
    ${TT} samtools view -hbo ${PFX}-realigned.bam
fi

# Merge and sort
if [ ! -s ${PFX}-final.bam ]; then
    ${TT} samtools cat ${PFX}-committed.bam ${PFX}-realigned.bam | \
    ${TT} samtools sort -@ ${THR} -o ${PFX}-final.bam
fi

