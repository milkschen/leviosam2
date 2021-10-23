set -x

ALN_RG=""
THR=$(nproc)
LEVIOSAM=leviosam
TIME=time # GNU time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step
FRAC_CLIPPED=0.95
ISIZE=1000
DEFERRED_BED=""
DISCARD_BED=""

# Default parameters for Bowtie 2
# ALN=bowtie2
# MAPQ=10
# ALN_SCORE=-10

# Default parameters for BWA-MEM
ALN=bwamem
MAPQ=30
ALN_SCORE=100


while getopts a:A:b:C:D:F:i:L:M:o:q:r:t:T: flag
do
    case "${flag}" in
        a) ALN=${OPTARG};;
        A) ALN_SCORE=${OPTARG};;
        b) ALN_IDX=${OPTARG};;
        C) CLFT=${OPTARG};;
        D) DISCARD_BED=${OPTARG};;
        F) DEFERRED_BED=${OPTARG};;
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
echo "BED where reads get deferred: ${DEFERRED_BED}";
echo "BED where reads get discarded: ${DISCARD_BED}";
TH=$(( ${THR} * 2/3 ))
STH=$(( ${THR} - ${TH} ))
echo "Num. threads: ${THR} (${TH}/${STH})";

if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit
fi

if [ ! -f ${PFX}-committed.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log \
            ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
            -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED}
    else
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam \
        -S mapq,isize,aln_score,clipped_frac -M ${MAPQ} -A ${ALN_SCORE} -Z ${ISIZE} -L ${FRAC_CLIPPED}
        # ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    fi
fi

if [ ! -f ${PFX}-committed-high_map.bam ]; then
    bedtools intersect -a ${PFX}-committed.bam -b ${DEFERRED_BED} -v > ${PFX}-committed-high_map.bam
    bedtools intersect -a ${PFX}-committed.bam -b ${DEFERRED_BED} > ${PFX}-committed-low_map.bam
fi

if [ ! -f ${PFX}-deferred-collate.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        samtools cat -o ${PFX}-deferred-merged.bam ${PFX}-deferred.bam ${PFX}-committed-low_map.bam
        ${TIME} -v -o collate.time_log \
            samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred-merged.bam
    else
        samtools cat ${PFX}-deferred-merged.bam ${PFX}-deferred.bam ${PFX}-committed-low_map.bam
        samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred-merged.bam
        # samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    fi
fi

if [ ! -f ${PFX}-deferred-R1.fq ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o to_fastq.time_log samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    else
        samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    fi
fi

if [ ! -f ${PFX}-deferred-re_aligned-paired.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${TIME} -v -o aln_paired.time_log bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        else
            ${TIME} -v -o aln_paired.time_log bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-R1.fq ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        fi
    else
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        else
            bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-R1.fq ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        fi
    fi
fi

if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${TIME} -v -o aln_unpaired.time_log bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        else
            ${TIME} -v -o aln_unpaired.time_log bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        fi
    else
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        else
            bwa mem -t ${THR} -R ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        fi
    fi
fi

if [ ! -f ${PFX}-merged.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge.time_log samtools cat -o ${PFX}-merged.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam ${PFX}-committed-high_map.bam
    else
        samtools cat -o ${PFX}-merged.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam ${PFX}-committed-high_map.bam
    fi
fi

if [ ! -f ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o sort_all.time_log samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    else
        samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    fi
fi
