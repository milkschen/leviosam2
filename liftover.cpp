#include <getopt.h>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include "liftover.hpp"

struct lift_opts {
    std::string vcf_fname = "";
    std::string sample = "";
    std::string outpre = "out";
    std::string lift_fname = "";
    std::string sam_fname = "";
};


lift::LiftMap lift_from_vcf(std::string fname, std::string sample);


void serialize_run(lift_opts args) {
    lift::LiftMap l(lift_from_vcf(args.vcf_fname, args.sample));
    std::ofstream o(args.outpre + ".lft");
    l.serialize(o);
    fprintf(stderr, "liftover file saved to %s\n", (args.outpre + ".lft").data());
}


void lift_run(lift_opts args) {
    // if "-l" not specified, then create a liftover
    lift::LiftMap l = [&]{
        if (args.lift_fname != "") {
            std::ifstream in(args.lift_fname);
            return lift::LiftMap(in);
        } else if (args.vcf_fname != "" && args.sample != "") {
            return lift::LiftMap(lift_from_vcf(args.vcf_fname, args.sample));
        } else {
            fprintf(stderr, "not enough parameters specified to build/load lift-over\n");
            exit(1);
        }
    } ();

    // read sam file here
    if (args.sam_fname == "") {
        fprintf(stderr, "no sam file specified! Use -a option\n");
        exit(1);
    }
    samFile* sam_fp = sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();

    // the "ref" lengths are all stored in the liftover structure. How do we loop over it?
    std::vector<std::string> contig_names;
    std::vector<size_t> contig_reflens;
    std::tie(contig_names, contig_reflens) = l.get_reflens();
    // for now we'll just write out the samfile raw
    // samFile* out_sam_fp = sam_open((args.outpre + ".sam").data(), "w");
    FILE* out_sam_fp = fopen((args.outpre + ".sam").data(), "w");
    fprintf(out_sam_fp, "@HD\tVN:1.6\tSO:unknown\n");
    for (auto i = 0; i < contig_names.size(); ++i) {
        fprintf(out_sam_fp, "@SQ\tSN:%s\tLN:%ld\n", contig_names[i].data(), contig_reflens[i]);
    }
    while (sam_read1(sam_fp, hdr, aln) >= 0) {
        bam1_core_t c = aln->core;
        std::string ref_name(hdr->target_name[c.tid]);
        fprintf(out_sam_fp, "%s\t", bam_get_qname(aln));
        fprintf(out_sam_fp, "%d\t", c.flag);
        if (c.flag & 4) { // unmapped here
            fprintf(out_sam_fp, "*\t");
            fprintf(out_sam_fp, "*\t");
            fprintf(out_sam_fp, "255\t"); // set MAPQ to unknown (255)
            fprintf(out_sam_fp, "*\t");
            fprintf(out_sam_fp, "*\t"); // RNEXT
            fprintf(out_sam_fp, "0\t"); // PNEXT
            fprintf(out_sam_fp, "0\t"); // TLEN (can probably copy?)
        } else {
            fprintf(out_sam_fp, "%s\t", ref_name.data()); // REF NAME
            /**** LIFTOVER STEP ****/
            fprintf(out_sam_fp, "%ld\t", l.alt_to_ref(ref_name, c.pos) + 1);  // POS
            /****               ****/
            fprintf(out_sam_fp, "255\t"); // set MAPQ to unknown (255)
            fprintf(out_sam_fp, "%s\t", l.cigar_alt_to_ref(ref_name, aln).data()); // CIGAR
            fprintf(out_sam_fp, "*\t"); // RNEXT
            fprintf(out_sam_fp, "0\t"); // PNEXT
            fprintf(out_sam_fp, "0\t"); // TLEN (can probably copy?)
        }
        // get query sequence
        std::string query_seq("");
        uint8_t* seq = bam_get_seq(aln);
        for (auto i = 0; i < c.l_qseq; ++i) {
            query_seq += seq_nt16_str[bam_seqi(seq, i)];
        }
        fprintf(out_sam_fp, "%s\t", query_seq.data());
        // get quality
        std::string qual_seq("");
        uint8_t* qual = bam_get_qual(aln);
        if (qual[0] == 255) qual_seq = "*";
        else {
            for (auto i = 0; i < c.l_qseq; ++i) {
                qual_seq += (char) (qual[i] + 33);
            }
        }
        fprintf(out_sam_fp, "%s", qual_seq.data());
        // TODO reconcile any tags that also need to be added.
        fprintf(out_sam_fp, "\n");
    }
    fclose(out_sam_fp);
}


lift::LiftMap lift_from_vcf(std::string fname, std::string sample) {
    if (fname == "" && sample == "") {
        fprintf(stderr, "vcf file name and sample name are required!! \n");
        exit(1);
    }
    vcfFile* fp = bcf_open(fname.data(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    return lift::LiftMap(fp, hdr, sample);
}


int main(int argc, char** argv) {
    int c;
    lift_opts args;
    static struct option long_options[] {
        {"vcf", required_argument, 0, 'v'},
        {"sample", required_argument, 0, 's'},
        {"prefix", required_argument, 0, 'p'},
        {"liftover", required_argument, 0, 'l'},
        {"sam", required_argument, 0, 'a'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "v:s:p:l:a:", long_options, &long_index)) != -1) { 
        switch (c) {
            case 'v':
                args.vcf_fname = optarg;
                break;
            case 's':
                args.sample = optarg;
                break;
            case 'p':
                args.outpre = optarg;
                break;
            case 'l':
                args.lift_fname = optarg;
                break;
            case 'a':
                args.sam_fname = optarg;
                break;
            default:
                fprintf(stderr, "ignoring option %s\n", c);
                exit(1);
                break;
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }

    if (!strcmp(argv[optind], "lift")) {
        lift_run(args);
    } else if (!strcmp(argv[optind], "serialize")) {
        serialize_run(args);
    }
    return 0;
}
