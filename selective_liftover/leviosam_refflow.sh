set -x

ALN_RG=""
THR=$(nproc)
MAPQ=10
LEVIOSAM=leviosam
TIME=time # GNU time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step
ALN=bowtie2

while getopts i:o:b:r:C:t:q:a:L:T:M: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        o) PFX=${OPTARG};;
        b) ALN_IDX=${OPTARG};;
        r) ALN_RG=${OPTARG};;
        C) CLFT=${OPTARG};;
        t) THR=${OPTARG};;
        q) MAPQ=${OPTARG};;
        a) ALN=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        T) TIME=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Aligner: ${ALN}";
echo "Aligner indexes prefix: ${ALN_IDX}";
echo "Aligner read group: ${ALN_RG}";
echo "LevioSAM index: ${CLFT}";
echo "LevioSAM MAPQ cutoff: ${MAPQ}";
TH=$(( ${THR} * 2/3 ))
STH=$(( ${THR} - ${TH} ))
echo "Num. threads: ${THR} (${TH}/${STH})";

if [[ ! ${ALN} =~ ^(bowtie2|bwamem)$ ]]; then
    echo "Invalid ${ALN}. Accepted input: bowtie2, bwamem"
    exit
fi

if [ ! -f ${PFX}.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    else
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    fi
fi

if [ ! -f ${PFX}-deferred-collate.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o collate.time_log samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    else
        samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
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
            ${TIME} -v -o aln_paired.time_log bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-R1.fq ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        fi
    else
        # bowtie2 ${ALN_RG} -p ${TH} -x ${ALN_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        else
            bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-R1.fq ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
        fi
    fi
fi
if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        if [[ ${ALN} == "bowtie2" ]]; then
            ${TIME} -v -o aln_unpaired.time_log bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        else
            ${TIME} -v -o aln_unpaired.time_log bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        fi
    else
        # bowtie2 ${ALN_RG} -p ${TH} -x ${ALN_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
        if [[ ${ALN} == "bowtie2" ]]; then
            bowtie2 ${ALN_RG} -p ${THR} -x ${ALN_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        else
            bwa mem -t ${THR} ${ALN_RG} ${ALN_IDX} ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        fi
    fi
fi

if [ ! -f ${PFX}-merged.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o merge.time_log samtools cat -o ${PFX}-merged.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam ${PFX}.bam
    else
        samtools cat -o ${PFX}-merged.bam ${PFX}-deferred-re_aligned-paired.bam ${PFX}-deferred-re_aligned-unpaired.bam ${PFX}.bam
    fi
fi

if [ ! -f ${PFX}-final.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o sort_all.time_log samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    else
        samtools sort -@ ${THR} -o ${PFX}-final.bam ${PFX}-merged.bam
    fi
fi
