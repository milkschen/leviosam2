# leviosam2.sh
#
# Full levioSAM2 pipeline for a SAM/BAM file
#
# Author: Nae-Chyun Chen
# Dept. of Computer Science, Johns Hopkins University
#
# Distributed under the MIT license
# https://github.com/milkschen/leviosam2
#

# Path to levioSAM2
LEVIOSAM=leviosam2
# Path to GNU time (use along with `MEASURE_TIME`)
TIME=time
# Set to 1 to measure time for each step
MEASURE_TIME=0
# Set to 1 to keep tmp files
KEEP_TMP=0
# Number of threads to use
THR=$(nproc)
# Set to 1 to use single-end mode
SINGLE_END=0

SOURCE_LABEL="source"
TARGET_LABEL="target"
ALN_RG=""

# LevioSAM2 parameters
ALLOWED_GAPS=0
FRAC_CLIPPED="-S clipped_frac:0.95"
ISIZE="-S isize:1000"
HDIST="-S hdist:5"
DEFER_DEST_BED=""
COMMIT_SOURCE_BED=""
REALN_CONFIG=""

# Default parameters for Bowtie 2 (local)
ALN=bowtie2
MAPQ="-S mapq:10"
ALN_SCORE="-S aln_score:-10"

# Default parameters for BWA-MEM
# ALN=bwamem
# MAPQ="-S mapq:30"
# ALN_SCORE="-S aln_score:100"

function usage {
    echo "Run full levioSAM2 pipeline for a SAM/BAM file"
    echo "Usage: $0 [-h] [options] -i <sam/bam> -o <prefix> -b <target.idx> -C <clft> -f <target.fasta>"
    echo "Author: Nae-Chyun Chen"
    echo ""
    echo "LevioSAM2 is distributed under the MIT license (2022)"
    echo "GitHub: https://github.com/milkschen/leviosam2"
    echo ""
    echo "Inputs:"
    echo "  -i path     Path to the input SAM/BAM file"
    echo "  -o prefix   Prefix of output files"
    echo "  -b string   Prefix to the aligner index"
    echo "  -C path     Path to the levioSAM2 index"
    echo "  -f path     Path to the target FASTA file"
    echo ""
    echo "Options:"
    echo "    -h          Display help message"
    echo "  System & Pipeline:"
    echo "    -K          Toggle to keep tmp files [off]"
    echo "    -M          Toggle to measure time using GNU time [off]"
    echo "    -S          Toggle to use single-end mode [off]"
    echo "    -t INT      Number of threads to use [max]"
    echo "    -L path     Path to the levioSAM2 binary [leviosam2]"
    echo "    -T path     Path to the GNU time binary [time]"
    echo "  LevioSAM2-lift:"
    echo "    -g INT      Number of gaps allowed during leviosam2-lift [0]"
    echo "    -X path     Path to the levioSAM2 re-alignment config YAML []"
    echo "  Commit/Defer/Suppress rules:"
    echo "    -A INT      Alignment score cutoff for the defer rule [-10]"
    echo "    -H INT      Edit distance cutoff for the defer rule [5]"
    echo "    -q INT      MAPQ cutoff for the defer rule [10]"
    echo "    -D path     Path to the force-defer annotation [empty]"
    echo "  Aligner:"
    echo "    -a string   Aligner to use (bowtie2|bwamem|minimap2|winnowmap2) [bowtie2]"
    echo "    -r string   The read group (RG) string []"
    echo "    -R path     Path to the suppress annotation []"
    exit 0
}

while getopts hKSPa:A:b:C:D:f:H:g:i:L:o:q:r:R:t:T:x: flag
do
    case "${flag}" in
        h) usage;;
        K) KEEP_TMP=1;;
        M) MEASURE_TIME=1;;
        S) SINGLE_END=1;;
        a) ALN=${OPTARG};;
        A) ALN_SCORE=" -S aln_score:${OPTARG}";;
        b) ALN_IDX=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFER_DEST_BED=" -D ${OPTARG}";;
        f) REF=${OPTARG};;
        H) HDIST=" -S hdist:${OPTARG}";;
        g) ALLOWED_GAPS=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        o) PFX=${OPTARG};;
        q) MAPQ=" -S mapq:${OPTARG}";;
        r) ALN_RG=${OPTARG};;
        R) COMMIT_SOURCE_BED=" -r ${OPTARG}";;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
        x) REALN_CONFIG=" -x ${OPTARG}";;
    esac
done

if [[ ${INPUT} == "" ]]; then
    echo "Input is not set"
    exit 1
fi
if [[ ${PFX} == "" ]]; then
    echo "Prefix is not set"
    exit 1
fi
if [[ ${REF} == "" ]]; then
    echo "Targer reference is not set"
    exit 1
fi
if [[ ${ALN_IDX} == "" ]]; then
    echo "Targer reference index is not set"
    exit 1
fi

MT=""
if (( ${MEASURE_TIME} > 0 )); then
    MT="${TIME} -v -ao leviosam2.time_log "
fi

if [[ ! ${ALN} =~ ^(bowtie2|bwamem|minimap2|winnowmap2)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem, minimap2, winnowmap2"
    exit 1
fi

set -xp

# Lifting over using leviosam2
if [ ! -s ${PFX}-committed.bam ]; then
    ${MT} ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
    ${MAPQ} ${ISIZE} ${ALN_SCORE} ${FRAC_CLIPPED} ${HDIST} \
    -G ${ALLOWED_GAPS} \
    ${REALN_CONFIG} \
    -m -f ${REF} ${DEFER_DEST_BED} ${COMMIT_SOURCE_BED}
fi

if (( ${SINGLE_END} == 1 )); then
    # Convert deferred reads to FASTQ
    if [ ! -s ${PFX}-deferred.fq.gz ]; then
        ${MT} samtools fastq ${PFX}-deferred.bam | \
        ${MT} bgzip > ${PFX}-deferred.fq.gz
    fi

    # Re-align deferred reads
    if [ ! -s ${PFX}-realigned.bam ]; then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${MT} bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} \
            -U ${PFX}-deferred.fq.gz |\
            ${MT} samtools view -hbo ${PFX}-realigned.bam
        elif [[ ${ALN} == "bwamem" ]]; then
            ALN_RG="-R ${ALN_RG}"
            ${MT} bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} \
            ${PFX}-deferred.fq.gz |\
            ${MT} samtools view -hbo ${PFX}-realigned.bam
        else
            ALN_RG="-R ${ALN_RG}"
            ${MT} ${ALN} -ax map-hifi --MD -t ${THR} ${ALN_RG} \
            ${REF} ${PFX}-deferred.fq.gz | \
            ${MT} samtools view -hbo ${PFX}-realigned.bam
        fi
    fi

    # Merge and sort
    if [ ! -s ${PFX}-final.bam ]; then
        ${MT} samtools cat ${PFX}-committed.bam ${PFX}-realigned.bam | \
        ${MT} samtools sort -@ ${THR} -o ${PFX}-final.bam
    fi

    # Clean tmp files
    if (( ${KEEP_TMP} < 1 )); then
        rm ${PFX}-committed.bam ${PFX}-deferred.bam 
    fi
else
    # Collate
    if [ ! -s ${PFX}-paired-deferred-R1.fq.gz ]; then
        ${MT} ${LEVIOSAM} collate \
        -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
    fi

    # Re-align deferred reads
    if [ ! -s ${PFX}-paired-realigned.bam ]; then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${MT} bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} \
            -1 ${PFX}-paired-deferred-R1.fq.gz -2 ${PFX}-paired-deferred-R2.fq.gz | \
            ${MT} samtools view -hb > ${PFX}-paired-realigned.bam
        elif [[ ${ALN} == "bwamem" ]]; then
            ALN_RG="-R ${ALN_RG}"
            ${MT} bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} \
            ${PFX}-paired-deferred-R1.fq.gz ${PFX}-paired-deferred-R2.fq.gz | \
            ${MT} samtools view -hb > ${PFX}-paired-realigned.bam
        else
            echo "We do not support paired-end mode for aligner ${ALN}"
            exit 1
        fi
    fi

    # Reference flow-style merging
    if [ ! -s ${PFX}-paired-deferred-reconciled.bam ]; then
        ${MT} samtools sort -@ ${THR} -n \
            -o ${PFX}-paired-realigned-sorted_n.bam ${PFX}-paired-realigned.bam
        ${MT} samtools sort -@ ${THR} -n \
            -o ${PFX}-paired-deferred-sorted_n.bam ${PFX}-paired-deferred.bam
        ${MT} ${LEVIOSAM} reconcile \
            -s ${SOURCE_LABEL}:${PFX}-paired-deferred-sorted_n.bam \
            -s ${TARGET_LABEL}:${PFX}-paired-realigned-sorted_n.bam \
            -m -o ${PFX}-paired-deferred-reconciled.bam
    fi

    # Merge, sort, and clean
    if [ ! -s ${PFX}-final.bam ]; then
        ${MT} samtools cat ${PFX}-paired-committed.bam ${PFX}-paired-deferred-reconciled.bam | \
        ${MT} samtools sort -@ ${THR} -o ${PFX}-final.bam
    fi

    # Clean tmp files
    if (( ${KEEP_TMP} < 1 )); then
        rm ${PFX}-paired-deferred.bam ${PFX}-paired-deferred-sorted_n.bam 
        rm ${PFX}-paired-realigned.bam ${PFX}-paired-realigned-sorted_n.bam
        rm ${PFX}-paired-deferred-reconciled.bam
        rm ${PFX}-paired-committed.bam ${PFX}-deferred.bam
    fi
fi
