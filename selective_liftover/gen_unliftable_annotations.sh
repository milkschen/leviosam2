set -xp 

DIR_LEVIOSAM="/home/cnaechy1/data_blangme2/naechyun/levioSAM" #$(pwd)
# Inputs

# GRCh38
# CHAIN="/home/cnaechy1/data_blangme2/fasta/chain/chm13_v1.1_hg2Y-grch38.chain"
# SOURCE_FA="/home/cnaechy1/data_blangme2/fasta/chm13_v1.1-hg2Y/chm13.v1.1-hg002_Y.fasta"
# DEST_FA="/home/cnaechy1/data_blangme2/fasta/grch38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
# # DEST_KMER_IDX needs to be consistent with `k`. Can leave this empty
# DEST_KMER_IDX=" -W /home/cnaechy1/data_blangme2/fasta/grch38/repetitive_k19.txt"
# PFX="chm13_v1.1_hg2Y-grch38"

# hg19
CHAIN="/home/cnaechy1/data_blangme2/fasta/chain/chm13_v1.1_hg2Y-hg19.chain"
SOURCE_FA="/home/cnaechy1/data_blangme2/fasta/chm13_v1.1-hg2Y/chm13.v1.1-hg002_Y.fasta"
DEST_FA="/home/cnaechy1/data_blangme2/fasta/grch37/hg19.fa"
# DEST_KMER_IDX needs to be consistent with `k`. Can leave this empty
DEST_KMER_IDX=""
PFX="chm13_v1.1_hg2Y-hg19"
MIN_SEG_SIZE="5000"

# Resulting files
SOURCE_BED="chm13.v1.1-hg002_Y.bed"

python ${DIR_LEVIOSAM}/scripts/verbosify_chain.py -c ${CHAIN} -o ${CHAIN}.annotated -b ${PFX}
python ${DIR_LEVIOSAM}/scripts/fai_to_bed.py -f ${SOURCE_FA}.fai -o ${SOURCE_BED}
bedtools sort -i ${PFX}-source.bed | \
    bedtools merge -i stdin | \
    bedtools subtract -a ${SOURCE_BED} -b stdin > ${PFX}-source-unliftable.bed
# bedtools sort -i ${PFX}-source.bed | bedtools merge -i stdin > ${PFX}-source-sorted_merged.bed
# bedtools subtract -a ${SOURCE_BED} -b ${PFX}-source-sorted_merged.bed > ${PFX}-source-unliftable.bed
python ${DIR_LEVIOSAM}/scripts/filter_bed_by_size.py -b ${PFX}-source-unliftable.bed \
    -o ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}.bed -s ${MIN_SEG_SIZE}


python ${DIR_LEVIOSAM}/scripts/extract_seq_from_bed.py -f ${SOURCE_FA} \
    -b ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}.bed \
    -o ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}.fasta
winnowmap -t $(nproc) -ak 19 ${DEST_KMER_IDX} ${DEST_FA} \
    ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}.fasta > \
    ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}-winnowmap2.sam

python ${DIR_LEVIOSAM}/scripts/sam_qname_to_bed.py -s ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}-winnowmap2.sam \
    -m 1-5 -o ${PFX}-source-unliftable-s_${MIN_SEG_SIZE}-winnowmap2-exc_1_5.bed
