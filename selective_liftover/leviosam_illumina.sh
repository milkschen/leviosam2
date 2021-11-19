set -xp

ALN_RG=""
THR=$(nproc)
LEVIOSAM=leviosam
TIME=time # GNU time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step
FRAC_CLIPPED=0.95
ISIZE=1000
DEFERRED_DEST_BED=""
DISCARD_BED=""

# Default parameters for Bowtie 2 (local)
# ALN=bowtie2
# MAPQ=10
# ALN_SCORE=200

# Default parameters for BWA-MEM
ALN=bwamem
MAPQ=30
ALN_SCORE=100


while getopts a:A:b:C:D:i:L:M:o:q:r:R:t:T: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        A) ALN_SCORE=${OPTARG};;
        b) ALN_IDX=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DEFERRED_DEST_BED=${OPTARG};;
        R) DISCARD_BED=${OPTARG};;
        i) INPUT=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
        o) PFX=${OPTARG};;
        q) MAPQ=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        t) THR=${OPTARG};;
        T) TIME=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Aligner: ${ALN}";
echo "Aligner indexes prefix: ${ALN_IDX}";
echo "Aligner read group: ${ALN_RG}";
echo "LevioSAM software: ${LEVIOSAM}";
echo "LevioSAM index: ${CLFT}";
echo "LevioSAM min MAPQ: ${MAPQ}";
echo "LevioSAM min AS: ${ALN_SCORE}";
echo "BED where reads get deferred: ${DEFERRED_DEST_BED}";
echo "BED where reads get discarded: ${DISCARD_BED}";
TH=$(( ${THR} * 2/3 ))
STH=$(( ${THR} - ${TH} ))
echo "Num. threads: ${THR} (${TH}/${STH})";

if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit
fi

# Lifting over using leviosam
if [ ! -f ${PFX}-committed.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${DEFERRED_DEST_BED} == "" ]]; then
            ${TIME} -v -o lift.time_log \
                ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
                -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0
        else
            ${TIME} -v -o lift.time_log \
                ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
                -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0 \
                -D ${DEFERRED_DEST_BED}
        fi
    else
        if [[ ${DEFERRED_DEST_BED} == "" ]]; then
            ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
            -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0
        else
            ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
            -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED} -G 0 \
            -D ${DEFERRED_DEST_BED}
        fi
    fi
fi

# Collate
if (( ${MEASURE_TIME} > 0)); then
    ${TIME} -v -o collate.time_log \
        ${LEVIOSAM} collate -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
else
    ${LEVIOSAM} collate -a ${PFX}-committed.bam -b ${PFX}-deferred.bam -p ${PFX}-paired
fi

# Re-align deferred reads
if [ ! -f ${PFX}-paired-realigned.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${TIME} -v -o aln_paired.time_log bowtie2 --local ${ALN_RG} -p ${THR} -x ${ALN_IDX} -1 ${PFX}-paired-deferred-R1.fq -2 ${PFX}-paired-deferred-R2.fq | samtools view -hb > ${PFX}-paired-realigned.bam
        else
            ${TIME} -v -o aln_paired.time_log bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-paired-deferred-R1.fq ${PFX}-paired-deferred-R2.fq | samtools view -hb > ${PFX}-paired-realigned.bam
        fi
    else
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 --local ${ALN_RG} -p ${THR} -x ${ALN_IDX} -1 ${PFX}-paired-deferred-R1.fq -2 ${PFX}-paired-deferred-R2.fq | samtools view -hb > ${PFX}-paired-realigned.bam
        else
            bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-paired-deferred-R1.fq ${PFX}-paired-deferred-R2.fq | samtools view -hb > ${PFX}-paired-realigned.bam
        fi
    fi
fi

# Merge and sort
if [ ! -f ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge.time_log \
            samtools cat -o ${PFX}-merged.bam ${PFX}-paired-committed.bam ${PFX}-paired-realigned.bam
        ${TIME} -v -o sort_all.time_log \
            samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    else
        samtools cat -o ${PFX}-merged.bam ${PFX}-paired-committed.bam ${PFX}-paired-realigned.bam
        samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    fi
fi
