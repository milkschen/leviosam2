#include <cstdio>
#include <map>
#include <stdio.h>
#include <thread>
#include <tuple>
#include <vector>
#include <zlib.h>
#include <getopt.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include "leviosam.hpp"

KSEQ_INIT(gzFile, gzread)
;;

extern "C" {
    int bam_fillmd1(bam1_t *b, const char *ref, int flag, int quiet_mode);
}

struct lift_opts {
    std::string vcf_fname = "";
    std::string chain_fname = "";
    std::string sample = "";
    std::string outpre = "";
    std::string out_format = "sam";
    std::string lift_fname = "";
    std::string sam_fname = "";
    std::string cmd = "";
    std::string haplotype = "0";
    int threads = 1;
    int chunk_size = 256;
    int verbose = 0;
    NameMap name_map;
    LengthMap length_map;
    int md_flag = 0;
    std::string ref_name = "";
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

LengthMap parse_length_map(const char* fname) {
    LengthMap lengths;
    FILE* fp = fopen(fname, "r");
    int x;
    char n[255];
    char l[255];
    while ((x = fscanf(fp, "%[^\t]\t%[^\n]\n", n, l)) != EOF) {
        if (x <= 0 || lengths.find(n) != lengths.end()) {
            fprintf(stderr, "error encountered reading length map\n");
            exit(1);
        }
        lengths[n] = std::atoi(l);
    }
    fclose(fp);
    return lengths;
}

void serialize_run(lift_opts args) {
    lift::LiftMap l(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map, args.length_map));
    if (args.outpre == "") {
        fprintf(stderr, "no output prefix specified! Use -p \n");
        print_serialize_help_msg();
        exit(1);
    }
    std::ofstream o(args.outpre + ".lft", std::ios::binary);
    l.serialize(o);
    fprintf(stderr, "levioSAM file saved to %s\n", (args.outpre + ".lft").data());
}

void read_and_lift(
    std::mutex* mutex_fread,
    std::mutex* mutex_fwrite,
    std::mutex* mutex_vec,
    samFile* sam_fp,
    samFile* out_sam_fp,
    bam_hdr_t* hdr,
    int chunk_size,
    lift::LiftMap* l,
    std::vector<std::string>* chroms_not_found,
    std::map<std::string, std::string>* ref_dict,
    int md_flag
){
    std::string ref_name;
    std::string ref;
    std::vector<bam1_t*> aln_vec;
    for (int i = 0; i < chunk_size; i++){
        bam1_t* aln = bam_init1();
        aln_vec.push_back(aln);
    }
    int read = 1;
    while (read >= 0){
        int num_actual_reads = chunk_size;
        {
            // read from SAM, protected by mutex
            std::lock_guard<std::mutex> g(*mutex_fread);
            for (int i = 0; i < chunk_size; i++){
                read = sam_read1(sam_fp, hdr, aln_vec[i]);
                if (read < 0){
                    num_actual_reads = i;
                    break;
                }
            }
        }
        for (int i = 0; i < num_actual_reads; i++){
            bam1_t* aln = aln_vec[i];
            bam1_core_t c = aln->core;
            std::string s1_name, s2_name;
            size_t pos;
            // If a read is mapped, lift its position.
            // If a read is unmapped, but its mate is mapped, lift its mates position.
            // Otherwise leave it unchanged.
            if (!(c.flag & BAM_FUNMAP)) {
                s2_name = hdr->target_name[c.tid];
                s1_name = l->get_other_name(s2_name);
                // chroms_not_found needs to be protected because chroms_not_found is shared
                pos = l->s2_to_s1(s2_name, c.pos, chroms_not_found, mutex_vec);
                // CIGAR
                l->cigar_s2_to_s1(s2_name, aln);
                aln->core.pos = pos;
            } else if ((c.flag & BAM_FPAIRED) && !(c.flag & BAM_FMUNMAP)) {
                s2_name = hdr->target_name[c.mtid];
                s1_name = l->get_other_name(s2_name);
                // chroms_not_found needs to be protected because chroms_not_found is shared
                pos = l->s2_to_s1(s2_name, c.mpos, chroms_not_found, mutex_vec);
                aln->core.pos = pos;
            }
            // Lift mate position.
            if (c.flag & BAM_FPAIRED) {
                // If the mate is unmapped, use the lifted position of the read.
                // If the mate is mapped, lift its position.
                if (c.flag & BAM_FMUNMAP){
                    aln->core.mpos = aln->core.pos;
                } else {
                    std::string ms2_name(hdr->target_name[c.mtid]);
                    std::string ms1_name(l->get_other_name(ms2_name));
                    size_t mpos = l->s2_to_s1(ms2_name, c.mpos, chroms_not_found, mutex_vec);
                    aln->core.mpos = mpos;
                    int isize = (c.isize == 0)? 0 : c.isize + (mpos - c.mpos) - (pos - c.pos);
                    aln->core.isize = isize;
                }
            }
            if (md_flag) {
                // change ref if needed
                if (s1_name != ref_name) {
                    ref = (*ref_dict)[s1_name];
                }
                bam_fillmd1(aln, ref.data(), md_flag, 1);
            }
            else { // strip MD and NM tags if md_flag not set bc the liftover invalidates them
                uint8_t* ptr = NULL;
                if ((ptr = bam_aux_get(aln, "MD")) != NULL) {
                    bam_aux_del(aln, ptr);
                }
                if ((ptr = bam_aux_get(aln, "NM")) != NULL) {
                    bam_aux_del(aln, bam_aux_get(aln, "NM"));
                }
            }
            ref_name = s1_name;
        }
        {
            // write to file, thread corruption protected by lock_guard
            std::lock_guard<std::mutex> g(*mutex_fwrite);
            // std::thread::id this_id = std::this_thread::get_id();
            for (int i = 0; i < num_actual_reads; i++){
                auto flag_write = sam_write1(out_sam_fp, hdr, aln_vec[i]);
                if (flag_write < 0){
                    std::cerr << "[Error] Failed to write record " << bam_get_qname(aln_vec[i]) << "\n";
                    exit(1);
                }
            }
        }
    }
    for (int i = 0; i < chunk_size; i++){
        bam_destroy1(aln_vec[i]);
    }
    aln_vec.clear();
}

std::map<std::string, std::string> load_fasta(std::string ref_name) {
    std::cerr << "Loading FASTA...";
    std::map<std::string, std::string> fmap;
    gzFile fp = gzopen(ref_name.data(), "r");
    kseq_t *seq;
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        fmap[std::string(seq->name.s)] = std::string(seq->seq.s);
    }
    kseq_destroy(seq);
    gzclose(fp);
    std::cerr << "done\n";
    return fmap;
}

void lift_run(lift_opts args) {
    std::cerr << "Loading levioSAM index...";
    // if "-l" not specified, then create a levioSAM
    lift::LiftMap l = [&]{
        if (args.lift_fname != "") {
            std::ifstream in(args.lift_fname, std::ios::binary);
            return lift::LiftMap(in);
        } else if (args.vcf_fname != "") {
            return lift::LiftMap(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map, args.length_map));
        } else {
            fprintf(stderr, "Not enough parameters specified to build/load lift-over\n");
            print_lift_help_msg();
            exit(1);
        }
    } ();

    std::cerr << "done\n";

    samFile* sam_fp = (args.sam_fname == "")?
        sam_open("-", "r") : sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* hdr = sam_hdr_read(sam_fp);

    // the "ref" lengths are all stored in the levio structure. How do we loop over it?
    std::vector<std::string> contig_names;
    std::vector<size_t> contig_reflens;
    std::tie(contig_names, contig_reflens) = l.get_s1_lens();

    std::string out_mode = (args.out_format == "sam")? "w" : "wb";
    samFile* out_sam_fp = (args.outpre == "-" || args.outpre == "")?
        sam_open("-", out_mode.data()) :
        sam_open((args.outpre + "." + args.out_format).data(), out_mode.data());

    // Write SAM headers. We update contig lengths if needed.
    sam_hdr_add_pg(hdr, "leviosam", "VN", VERSION, "CL", args.cmd.data(), NULL);
    for (auto i = 0; i < hdr->n_targets; i++) {
        auto contig_itr = std::find(contig_names.begin(), contig_names.end(), hdr->target_name[i]);
        // If a contig is in contig_names, look up the associated length.
        if (contig_itr != contig_names.end()) {
            auto len_str = std::to_string(contig_reflens[contig_itr - contig_names.begin()]);
            if (sam_hdr_update_line(
                    hdr, "SQ",
                    "SN", contig_names[contig_itr - contig_names.begin()].data(),
                    "LN", len_str.data(), NULL) < 0) {
                std::cerr << "[Error] failed when converting contig length for "
                          << hdr->target_name[i] << "\n";
                exit(1);
            }
        }
    }
    auto write_hdr = sam_hdr_write(out_sam_fp, hdr);


    std::map<std::string, std::string> ref_dict;
    if (args.md_flag) {
        // load
        if (args.ref_name == "") {
            std::cerr << "error: -m/--md -f <fasta> to be provided as well\n";
            exit(1);
        }
        ref_dict = load_fasta(args.ref_name);
    }

    // Store chromosomes found in SAM but not in the VCF.
    // We use a vector to avoid printing out the same warning msg multiple times.
    std::vector<std::string> chroms_not_found;
    const int num_threads = args.threads;
    const int chunk_size = args.chunk_size;
    std::vector<std::thread> threads;
    std::mutex mutex_fread, mutex_fwrite, mutex_vec;
    for (int j = 0; j < num_threads; j++){
        threads.push_back(
            std::thread(
                read_and_lift,
                &mutex_fread,
                &mutex_fwrite,
                &mutex_vec,
                sam_fp,
                out_sam_fp,
                hdr,
                chunk_size,
                &l,
                &chroms_not_found,
                &ref_dict,
                args.md_flag
            )
        );
    }
    for (int j = 0; j < num_threads; j++){
        if(threads[j].joinable())
            threads[j].join();
    }
    threads.clear();
    sam_close(out_sam_fp);
}

/*
lift::LiftMap lift_from_chain(std::string fname) {
    if (fname == "") {
        fprintf(stderr, "chain file name is required!! \n");
        print_serialize_help_msg();
        exit(1);
    }
    chain::ChainFile* cfp = chain::chain_open(fname);
    // chain::bcf_hdr_t* hdr = bcf_hdr_read(fp);
    return lift::LiftMap(cfp);
}
*/

lift::LiftMap lift_from_vcf(std::string fname, 
                            std::string sample, 
                            std::string haplotype, 
                            NameMap names, LengthMap lengths) {
    if (fname == "" && sample == "") {
        fprintf(stderr, "vcf file name and sample name are required!! \n");
        print_serialize_help_msg();
        exit(1);
    }
    vcfFile* fp = bcf_open(fname.data(), "r");
    bcf_hdr_t* hdr = bcf_hdr_read(fp);
    return lift::LiftMap(fp, hdr, sample, haplotype, names, lengths);
}


std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}


void print_serialize_help_msg(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   leviosam serialize [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "         -v string Build a leviosam index using a VCF file.\n");
    fprintf(stderr, "         -s string The sample used to build leviosam index (-v needs to be set).\n");
    fprintf(stderr, "         -p string The prefix of the output file.\n");
    fprintf(stderr, "         -O format Output file format, can be sam or bam. [sam]\n");
    fprintf(stderr, "         -g 0/1    The haplotype used to index leviosam. [0] \n");
    fprintf(stderr, "         -n string Path to a name map file.\n");
    fprintf(stderr, "                   This can be used to map '1' to 'chr1', or vice versa.\n");
    fprintf(stderr, "         -k string Path to a length map file.\n");
    fprintf(stderr, "\n");
}

void print_lift_help_msg(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Usage:   leviosam lift [options]\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "         -a string Path to the SAM/BAM file to be lifted. \n");
    fprintf(stderr, "                   Leave empty or set to \"-\" to read from stdin.\n");
    fprintf(stderr, "         -l string Path to a leviosam index.\n");
    fprintf(stderr, "         -v string If -l is not specified, can build indexes using a VCF file.\n");
    fprintf(stderr, "         -t INT    Number of threads used. [1] \n");
    fprintf(stderr, "         -T INT    Chunk size for each thread. [256] \n");
    fprintf(stderr, "                   Each thread queries <-T> reads, lifts, and writes.\n");
    fprintf(stderr, "                   Setting a higher <-T> uses slightly more memory but might benefit thread scaling.\n");
    fprintf(stderr, "         -m        add MD and NM to output alignment records (requires -f option)\n");
    fprintf(stderr, "         -f string Fasta reference that corresponds to input SAM/BAM (for use w/ -m option)\n");
    fprintf(stderr, "         The options for serialize can also be used here, if -v is set.\n");
    fprintf(stderr, "\n");
}

void print_main_help_msg(){
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: leviosam (lift over SAM/BAM based on VCF)\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Usage:   leviosam <command> [options]\n\n");
    fprintf(stderr, "Commands:serialize   Build a leviosam index.\n");
    fprintf(stderr, "         lift        Lift alignments.\n");
    fprintf(stderr, "Options: -h          Print detailed usage.\n");
    fprintf(stderr, "\n");
}

int main(int argc, char** argv) {
    int c;
    lift_opts args;
    args.cmd = make_cmd(argc,argv);
    static struct option long_options[] {
        {"vcf", required_argument, 0, 'v'},
        {"chain", required_argument, 0, 'c'},
        {"sample", required_argument, 0, 's'},
        {"prefix", required_argument, 0, 'p'},
        {"leviosam", required_argument, 0, 'l'},
        {"sam", required_argument, 0, 'a'},
        {"out_format", required_argument, 0, 'O'},
        {"haplotype", required_argument, 0, 'g'},
        {"threads", required_argument, 0, 't'},
        {"chunk_size", required_argument, 0, 'T'},
        {"md", required_argument, 0, 'm'},
        // {"nm", required_argument, 0, 'x'},
        {"reference", required_argument, 0, 'f'},
        {"verbose", no_argument, &args.verbose, 1},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "hmv:c:s:p:l:a:O:g:n:k:t:T:f:", long_options, &long_index)) != -1) {
        switch (c) {
            case 'v':
                args.vcf_fname = optarg;
                break;
            case 'c':
                args.chain_fname = optarg;
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
            case 'O':
                args.out_format = optarg;
                break;
            case 'g':
                args.haplotype = optarg;
                break;
            case 'n':
                args.name_map = parse_name_map(optarg);
                break;
            case 'k':
                args.length_map = parse_length_map(optarg);
                break;
            case 't':
                args.threads = atoi(optarg);
                break;
            case 'T':
                args.chunk_size = atoi(optarg);
                break;
            case 'm':
                args.md_flag |= 8;
                args.md_flag |= 16;
                break;
            case 'f':
                args.ref_name = std::string(optarg);
                break;
            case 'h':
                print_serialize_help_msg();
                print_lift_help_msg();
                exit(1);
            default:
                fprintf(stderr, "ignoring option %c\n", c);
                exit(1);
        }
    }

    if (argc - optind < 1) {
        fprintf(stderr, "no argument provided\n");
        print_main_help_msg();
        exit(1);
    }
    if (args.haplotype != "0" && args.haplotype != "1"){
        fprintf(stderr, "invalid haplotype %s\n", args.haplotype.c_str());
        exit(1);
    }

    if (args.out_format != "sam" && args.out_format != "bam") {
        fprintf(stderr, "Not supported extension format %s\n", args.out_format.c_str());
        exit(1);
    }

    if (!strcmp(argv[optind], "lift")) {
        lift_run(args);
    } else if (!strcmp(argv[optind], "serialize")) {
        serialize_run(args);
    }
    return 0;
}

