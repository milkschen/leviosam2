#include "liftover.hpp"
#include <getopt.h>
#include <cstdio>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

using NameMap = std::vector<std::pair<std::string,std::string>>;

struct lift_opts {
    std::string vcf_fname = "";
    std::string sample = "";
    std::string outpre = "";
    std::string lift_fname = "";
    std::string sam_fname = "";
    std::string cmd = "";
    std::string haplotype = "0";
    NameMap name_map;
};

NameMap parse_name_map(const char* fname) {
    NameMap names;
    FILE* fp = fopen(fname, "r");
    char n1[255];
    char n2[255];
    int x;
    while ((x = fscanf(fp,"%[^\t]\t%[^\n]\n", n1, n2)) != EOF) {
        if (x <= 0) {
            fprintf(stderr, "error encountered reading name map\n");
            exit(1);
        }
        names.push_back(std::make_pair(std::string(n1), std::string(n2)));
    }
    fclose(fp);
    return names;
}

lift::LiftMap lift_from_vcf(std::string fname, 
                            std::string sample, 
                            std::string haplotype, 
                            NameMap names);

void serialize_run(lift_opts args) {
    lift::LiftMap l(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map));
    if (args.outpre == "") {
        fprintf(stderr, "no output prefix specified! Use -p \n");
        exit(1);
    }
    std::ofstream o(args.outpre + ".lft");
    l.serialize(o);
    fprintf(stderr, "liftover file saved to %s\n", (args.outpre + ".lft").data());
}


char* get_PG(bam_hdr_t* hdr) {
    char* hdr_txt = (char*) malloc(sizeof(char) * hdr->l_text + 1);
    char* buf = (char*) malloc(sizeof(char) * hdr->l_text);
    buf[0] = '\0';
    std::strcpy(hdr_txt, hdr->text);
    char* token = std::strtok(hdr_txt, "@");
    while (token != NULL) {
        if (token[0] == 'P' && token[1] == 'G') {
            strcat(buf, "@");
            strcat(buf, token);
        }
        token = std::strtok(NULL, "@");
    }
    free(hdr_txt);
    return buf;
}


void lift_run(lift_opts args) {
    // if "-l" not specified, then create a liftover
    lift::LiftMap l = [&]{
        if (args.lift_fname != "") {
            std::ifstream in(args.lift_fname);
            return lift::LiftMap(in);
        } else if (args.vcf_fname != "") {
            return lift::LiftMap(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map));
        } else {
            fprintf(stderr, "not enough parameters specified to build/load lift-over\n");
            exit(1);
        }
    } ();

    fprintf(stderr, "loaded liftmap!\n");

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
    std::tie(contig_names, contig_reflens) = l.get_s1_lens();
    // for now we'll just write out the samfile raw
    FILE* out_sam_fp;
    if (args.outpre == "") 
        out_sam_fp = stdout;
    else 
        out_sam_fp = fopen((args.outpre + ".sam").data(), "w");
    fprintf(out_sam_fp, "@HD\tVN:1.6\tSO:unknown\n");
    for (auto i = 0; i < contig_names.size(); ++i) {
        fprintf(out_sam_fp, "@SQ\tSN:%s\tLN:%ld\n", contig_names[i].data(), contig_reflens[i]);
    }
    /* This is only supported in the latest dev release of htslib as of July 22. 
     * add back in at next htslib release
    kstring_t s = KS_INITIALIZE;
    sam_hdr_find_line_id(hdr, "PG", NULL, NULL, &s);
    fprintf(out_sam_fp, "%s\n", s.s);
    ks_free(&s);
    */
    char* prev_pg = get_PG(hdr);
    fprintf(out_sam_fp, "%s", prev_pg);
    free(prev_pg);
    fprintf(out_sam_fp, "@PG\tID:liftover\tPN:liftover\tCL:\"%s\"\n", args.cmd.data());
    while (sam_read1(sam_fp, hdr, aln) >= 0) {
        bam1_core_t c = aln->core;
        fprintf(out_sam_fp, "%s\t", bam_get_qname(aln));
        fprintf(out_sam_fp, "%d\t", c.flag);
        if (c.flag & 4) { // unmapped here
            fprintf(out_sam_fp, "*\t"); // RNAME (String)
            fprintf(out_sam_fp, "0\t"); // POS (Int)
            fprintf(out_sam_fp, "255\t"); // set MAPQ to unknown (255)
            fprintf(out_sam_fp, "*\t");
            fprintf(out_sam_fp, "*\t"); // RNEXT
            fprintf(out_sam_fp, "0\t"); // PNEXT
            fprintf(out_sam_fp, "0\t"); // TLEN (can probably copy?)
        } else {
            std::string ref_name(hdr->target_name[c.tid]);
            fprintf(out_sam_fp, "%s\t", ref_name.data()); // REF NAME
            /**** LIFTOVER STEP ****/
            fprintf(out_sam_fp, "%ld\t", l.s2_to_s1(ref_name, c.pos) + 1);  // POS
            /****               ****/
            // fprintf(out_sam_fp, "255\t"); // set MAPQ to unknown (255)
            fprintf(out_sam_fp, "%d\t", c.qual);
            fprintf(out_sam_fp, "%s\t", l.cigar_s2_to_s1(ref_name, aln).data()); // CIGAR
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


lift::LiftMap lift_from_vcf(std::string fname, 
                            std::string sample, 
                            std::string haplotype, 
                            NameMap names) {
    if (fname == "" && sample == "") {
        fprintf(stderr, "vcf file name and sample name are required!! \n");
        exit(1);
    }
    vcfFile* fp = bcf_open(fname.data(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    return lift::LiftMap(fp, hdr, sample, haplotype, names);
}


std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}


int main(int argc, char** argv) {
    int c;
    lift_opts args;
    args.cmd = make_cmd(argc,argv);
    static struct option long_options[] {
        {"vcf", required_argument, 0, 'v'},
        {"sample", required_argument, 0, 's'},
        {"prefix", required_argument, 0, 'p'},
        {"liftover", required_argument, 0, 'l'},
        {"sam", required_argument, 0, 'a'},
        {"haplotype", required_argument, 0, 'g'}
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "v:s:p:l:a:g:n:", long_options, &long_index)) != -1) {
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
            case 'g':
                args.haplotype = optarg;
                break;
            case 'n':
                args.name_map = parse_name_map(optarg);
                break;
            default:
                fprintf(stderr, "ignoring option %c\n", c);
                exit(1);
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "no argument provided\n");
        exit(1);
    }
    if (args.haplotype != "0" && args.haplotype != "1"){
        fprintf(stderr, "invalid haplotype %s\n", args.haplotype);
        exit(1);
    }
    if (!strcmp(argv[optind], "lift")) {
        lift_run(args);
    } else if (!strcmp(argv[optind], "serialize")) {
        serialize_run(args);
    }
    return 0;
}
