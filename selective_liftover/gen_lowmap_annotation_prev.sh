CHAIN="/home/cnaechy1/data_blangme2/fasta/chain/hg19ToHg38.over.chain"
SOURCE_REF="/home/cnaechy1/data_blangme2/fasta/grch37/hg19.fa"
DEST_REF="/home/cnaechy1/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
PFX="hg19ToHg38"
SOURCE_HIGHMAP="hg19-k100_e0.01-high.bed"
SOURCE_LOWMAP="hg19-k100_e0.01-low.bed"
SOURCE_TO_DEST_LOWMAP="hg19_to_h38-k100_e0.01-low.bed"
DEST_HIGHMAP="grch38-k100_e0.01-high.bed"

DIR_LEVIOSAM="."
UCSC_LIFTOVER="/home/cnaechy1/data_blangme2/naechyun/software/ucsc_utilities/liftOver"

# Get low-identity chain segments
python ${DIR_LEVIOSAM}/scripts/verbosify_chain.py -c ${CHAIN} -o ${PFX}.verbose.chain \
    -b ${PFX} -s ${PFX}.summary -f1 ${SOURCE_REF} -f2 ${DEST_REF}
python ${DIR_LEVIOSAM}/scripts/get_low_identity_regions.py -s ${PFX}.summary \
    -l 0.98 -o ${PFX}.identity_lt98.hg38.bed


# Calculate mappability annotations
# genmap2 has not been released
# genmap2 -k 100 -e 0.01 -s 10 -bg ${SOURCE_REF} # hg19.fa -> hg19_100_0.01_10_1.0.bedgraph
# genmap2 -k 100 -e 0.01 -s 10 -bg ${DEST_REF}

if [ ! -f ${SOURCE_REF}.fai ]; then
    samtools faidx ${SOURCE_REF}
fi
python ${DIR_LEVIOSAM}/scripts/fai_to_bed.py -f ${SOURCE_REF}.fai -o ${SOURCE_REF}.fai.bed

SOURCE_BG="/home/cnaechy1/data_blangme2/fasta/grch37/hg19_100_0.01_10_1.0.bedgraph"
python ${DIR_LEVIOSAM}/scripts/get_mappable_regions.py -b ${SOURCE_BG} -k 100 -o ${SOURCE_HIGHMAP}
bedtools subtract -a ${SOURCE_REF}.fai.bed -b ${SOURCE_HIGHMAP} > ${SOURCE_LOWMAP}
${UCSC_LIFTOVER} ${SOURCE_LOWMAP} ${CHAIN} ${SOURCE_TO_DEST_LOWMAP} ${SOURCE_TO_DEST_LOWMAP}.unmapped

DEST_BG="/home/cnaechy1/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_100_0.01_10_1.0.bedgraph"
python ${DIR_LEVIOSAM}/scripts/get_mappable_regions.py -b ${DEST_BG} -k 100 -o ${DEST_HIGHMAP}

# Low-mappability in SOURCE; high-mappability in DEST
bedtools intersect -a ${SOURCE_TO_DEST_LOWMAP} -b ${DEST_HIGHMAP} > k100_e0.01-hg19_lowmap-grch38_highmap-grch38.bed


# Merge low-identity and mappability-reduced BED files
if [ ! -f ${DEST_REF}.fai ]; then
    samtools faidx ${DEST_REF}
fi
cat ${PFX}.identity_lt98.hg38.bed k100_e0.01-hg19_lowmap-grch38_highmap-grch38.bed |\
    bedtools intersect -a stdin -b ${DEST_REF}.fai |\
    bedtools sort -i stdin |\
    bedtools merge -i stdin > hg19-grch38-lt98-map_reduction_k100_e0.01-merged.bed
