function usage {
    echo "Generate a chain file from a VCF file."
    echo "usage: $0 [-h] [-f fasta] [-v vcf] [-p prefix]"
    echo "  -f fasta    Path to a FASTA file"
    echo "  -v vcf      Path to a gz-compressed and indexed VCF file"
    echo "  -p prefix   Prefix for output files"
    echo "  -h          Display help"
    exit 0
}


while getopts "hf:v:p:" OPTION; do
case ${OPTION} in
    h)  # Display Help
        usage
        ;;
    f)  # Set reference path
        REF=${OPTARG}
        ;;
    v)  # Set VCF path
        VCF=${OPTARG}
        ;;
    p)  # Set output prefix
        PREFIX=${OPTARG}
        ;;
esac
done

#REF=/net/langmead-bigmem-ib.bluecrab.cluster/storage2/naechyun/reference_flow/resources/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
#VCF=wg-maj.vcf.gz
#PREFIX=wg-maj
CONSENSUS=${PREFIX}.fa
RCHAIN=${PREFIX}-reverse.chain
CHAIN=${PREFIX}.chain
if [[ -f $RCHAIN ]]; then
    echo "${RCHAIN} already exists. Skip bcftools consensus."
else
    bcftools consensus -c $RCHAIN -f $REF -o $CONSENSUS $VCF
fi
awk '{ if ($0 ~ /^chain/) {print} else {print $1, $3, $2} }' $RCHAIN > $CHAIN
