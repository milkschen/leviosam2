# Generate genomic annotations with low-liftability for a pair-wise whole genome alignment
#
# Author: Nae-Chyun Chen (2021)
#

set -xp 

CHAIN="/home/cnaechy1/data_blangme2/fasta/chain/chm13_v1.1_hg2Y-grch38.chain"
SOURCE_REF="/home/cnaechy1/data_blangme2/fasta/chm13_v1.1-hg2Y/chm13.v1.1-hg002_Y.fasta"
DEST_REF="/home/cnaechy1/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
SOURCE_LABEL="chm13_v1.1_hg2Y"
DEST_LABEL="grch38"
PFX="${SOURCE_LABEL}-${DEST_LABEL}"

K="100"
E="0.01"
L="0.98"

SOURCE_BG="/home/cnaechy1/data_blangme2/fasta/chm13_v1.1-hg2Y/chm13.v1.1-hg002_Y_${K}_${E}_10_1.0.bedgraph"
DEST_BG="/home/cnaechy1/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set_100_0.01_10_1.0.bedgraph"

SOURCE_HIGHMAP="${SOURCE_LABEL}-k${K}_e${E}-high.bed"
SOURCE_LOWMAP="${SOURCE_LABEL}-k${K}_e${E}-low.bed"
SOURCE_TO_DEST_HIGHMAP="${SOURCE_LABEL}_to_${DEST_LABEL}-k${K}_e${E}-high.bed"
DEST_HIGHMAP="${DEST_LABEL}-k${K}_e${E}-high.bed"
DEST_LOWMAP="${DEST_LABEL}-k${K}_e${E}-low.bed"

DIR_LEVIOSAM="/home/cnaechy1/data_blangme2/naechyun/levioSAM/"
UCSC_LIFTOVER="/home/cnaechy1/data_blangme2/naechyun/software/ucsc_utilities/liftOver"

# Get low-identity chain segments
if [ ! -s ${PFX}.summary ]; then
    python ${DIR_LEVIOSAM}/scripts/verbosify_chain.py -c ${CHAIN} -o ${PFX}.verbose.chain \
        -b ${PFX} -s ${PFX}.summary -f1 ${SOURCE_REF} -f2 ${DEST_REF}
fi
if [ ! -s ${PFX}-identity_lt${L}-${DEST_LABEL}.bed ]; then
    python ${DIR_LEVIOSAM}/scripts/get_low_identity_regions.py -s ${PFX}.summary \
        -l ${L} -o ${PFX}-identity_lt${L}-${DEST_LABEL}.bed
fi


# Calculate mappability annotations
# genmap2 has not been released
# genmap2 -k ${K} -e ${E} -s 10 -bg ${SOURCE_REF} # hg19.fa -> hg19_${K}_${E}_10_1.0.bedgraph
# genmap2 -k ${K} -e ${E} -s 10 -bg ${DEST_REF}

if [ ! -s ${SOURCE_HIGHMAP} ]; then
    python ${DIR_LEVIOSAM}/scripts/get_mappable_regions.py -b ${SOURCE_BG} -k ${K} -o ${SOURCE_HIGHMAP}
fi
if [ ! -s ${SOURCE_TO_DEST_HIGHMAP} ]; then
    ${UCSC_LIFTOVER} ${SOURCE_HIGHMAP} ${CHAIN} ${SOURCE_TO_DEST_HIGHMAP} ${SOURCE_TO_DEST_HIGHMAP}.unmapped
fi

if [ ! -s ${DEST_HIGHMAP} ]; then
    python ${DIR_LEVIOSAM}/scripts/get_mappable_regions.py -b ${DEST_BG} -k ${K} -o ${DEST_HIGHMAP}
fi

if [ ! -s ${DEST_REF}.fai ]; then
    samtools faidx ${DEST_REF}
fi
if [ ! -s ${DEST_LOWMAP} ]; then
    python ${DIR_LEVIOSAM}/scripts/fai_to_bed.py -f ${DEST_REF}.fai -o ${DEST_REF}.fai.bed
    bedtools subtract -a ${DEST_REF}.fai.bed -b ${DEST_HIGHMAP} > ${DEST_LOWMAP}
fi

# Low-mappability in DEST; high-mappability in SOURCE
# This is wrt the DEST coordinates
if [ ! -s k${K}_e${E}-${SOURCE_LABEL}_highmap-${DEST_LABEL}_lowmap-${DEST_LABEL}.bed ]; then
    bedtools intersect -a ${SOURCE_TO_DEST_HIGHMAP} -b ${DEST_LOWMAP} > k${K}_e${E}-${SOURCE_LABEL}_highmap-${DEST_LABEL}_lowmap-${DEST_LABEL}.bed
fi


# Merge low-identity and mappability-reduced BED files
cat ${PFX}-identity_lt${L}-${DEST_LABEL}.bed k${K}_e${E}-${SOURCE_LABEL}_highmap-${DEST_LABEL}_lowmap-${DEST_LABEL}.bed |\
    bedtools intersect -a stdin -b ${DEST_REF}.fai.bed |\
    bedtools sort -i stdin |\
    bedtools merge -i stdin > ${SOURCE_LABEL}-${DEST_LABEL}-lt${L}-map_reduction_k${K}_e${E}.bed
