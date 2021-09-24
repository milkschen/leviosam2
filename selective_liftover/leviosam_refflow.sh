set -x

THR=$(nproc)
MAPQ=10
LEVIOSAM=leviosam
TIME=time # GNU time
MEASURE_TIME=1 # Set to a >0 value to measure time for each step

while getopts i:o:b:r:C:t:q:L:T:M: flag
do
    case "${flag}" in
        i) INPUT=${OPTARG};;
        o) PFX=${OPTARG};;
        b) BT2_IDX=${OPTARG};;
        r) BT2_RG=${OPTARG};;
        C) CLFT=${OPTARG};;
        t) THR=${OPTARG};;
        q) MAPQ=${OPTARG};;
        L) LEVIOSAM=${OPTARG};;
        T) TIME=${OPTARG};;
        M) MEASURE_TIME=${OPTARG};;
    esac
done

echo "Input BAM: ${INPUT}";
echo "Output prefix: ${PFX}";
echo "Bowtie2 indexes prefix: ${BT2_IDX}";
echo "Bowtie2 read group: ${BT2_RG}";
echo "LevioSAM index: ${CLFT}";
echo "LevioSAM MAPQ cutoff: ${MAPQ}";
TH=$(( ${THR} * 2/3 ))
STH=$(( ${THR} - ${TH} ))
echo "Num. threads: ${THR} (${TH}/${STH})";

if [ ! -f ${PFX}.bam ]; then
    lift () {
        ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    }
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o lift.time_log lift 
        #${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    else
        lift
        # ${LEVIOSAM} lift -C ${CLFT} -a ${INPUT} -t ${THR} -p ${PFX} -O bam -S mapq -M ${MAPQ}
    fi
fi

if [ ! -f ${PFX}-deferred-collate.bam ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o collate.time_log samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    else
        samtools collate -@ ${THR} -fo ${PFX}-deferred-collate.bam ${PFX}-deferred.bam
    fi
fi
# FQSIZE=$(stat -c%s "${PFX}-deferred-R1.fq")
# if (( ${FQSIZE} > 0 )); then
#     echo "Skip converting FASTQ"
# else
if [ ! -f ${PFX}-deferred-R1.fq ]; then
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o to_fastq.time_log samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    else
        samtools fastq -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq -s ${PFX}-deferred-S.fq ${PFX}-deferred-collate.bam
    fi
fi

if [ ! -f ${PFX}-deferred-re_aligned-paired.bam ]; then
    aln_paired () {
        bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
    }
    if (( ${MEASURE_TIME} > 0)); then
        ${TIME} -v -o aln_paired.time_log aln_paired
        # bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
    else
        aln_paired
        # bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-paired.bam
        # bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -1 ${PFX}-deferred-R1.fq -2 ${PFX}-deferred-R2.fq | samtools view -hb > ${PFX}-deferred-re_aligned-paired.bam
    fi
fi
if [ ! -f ${PFX}-deferred-re_aligned-unpaired.bam ]; then
    aln_unpaired () {
        bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
    }
    if (( ${MEASURE_TIME} > 0)); then
        # ${TIME} -v -o aln_unpaired.time_log bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
        ${TIME} -v -o aln_unpaired.time_log aln_unpaired
    else
        aln_unpaired
        # bowtie2 ${BT2_RG} -p ${TH} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools sort -@ ${STH} -o ${PFX}-deferred-re_aligned-unpaired.bam
        # bowtie2 ${BT2_RG} -p ${THR} -x ${BT2_IDX} -U ${PFX}-deferred-S.fq | samtools view -hb > ${PFX}-deferred-re_aligned-unpaired.bam
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
