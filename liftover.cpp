#include "liftover.hpp"
#include <vector>
#include <tuple>
#include <unordered_map>
#include <getopt.h>
#include <cstdio>
#include <thread>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/kstring.h>

using NameMap = std::vector<std::pair<std::string,std::string>>;
using LengthMap = std::unordered_map<std::string,size_t>;

struct lift_opts {
    std::string vcf_fname = "";
    std::string sample = "";
    std::string outpre = "";
    std::string lift_fname = "";
    std::string sam_fname = "";
    std::string cmd = "";
    std::string haplotype = "0";
    int threads = 1;
    int chunk_size = 64;
    int verbose = 0;
    NameMap name_map;
    LengthMap length_map;
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

lift::LiftMap lift_from_vcf(std::string fname, 
                            std::string sample, 
                            std::string haplotype, 
                            NameMap names, LengthMap lengths);

void serialize_run(lift_opts args) {
    lift::LiftMap l(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map, args.length_map));
    if (args.outpre == "") {
        fprintf(stderr, "no output prefix specified! Use -p \n");
        exit(1);
    }
    std::ofstream o(args.outpre + ".lft", std::ios::binary);
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

/* just gets AS for now till I get around to iterating through all the tags 
 */
std::string tag_to_string(const bam1_t* rec) {
    std::string tag_string("");
    uint8_t* aux = bam_aux_get(rec, "AS");
    if (aux != NULL) {
        tag_string += "AS:i:";
        tag_string += std::to_string(bam_aux2i(aux));
    }
    aux = bam_aux_get(rec, "NM");
    if (aux != NULL) {
        tag_string += "\tNM:i:";
        tag_string += std::to_string(bam_aux2i(aux));
    }
    return tag_string;
}

void read_and_lift(
    std::mutex* mutex,
    samFile* sam_fp,
    FILE* out_sam_fp,
    bam_hdr_t* hdr,
    int chunk_size,
    lift::LiftMap* l,
    std::vector<std::string>* chroms_not_found
){
    // std::cerr << std::this_thread::get_id() << "\n";
    std::vector<bam1_t*> aln_vec;
    for (int i = 0; i < chunk_size; i++){
        bam1_t* aln = bam_init1();
        aln_vec.push_back(aln);
    }
    int read = 1;
    while(read >= 0){
        int num_actual_reads = chunk_size;
        mutex->lock();
        for (int i = 0; i < chunk_size; i++){
            read = sam_read1(sam_fp, hdr, aln_vec[i]);
            if (read < 0){
                num_actual_reads = i;
                break;
            }
        }
        mutex->unlock();
        for (int i = 0; i < num_actual_reads; i++){
            std::string sam_out;
            bam1_t* aln = aln_vec[i];
            bam1_core_t c = aln->core;
            sam_out += bam_get_qname(aln);
            sam_out += "\t";
            sam_out += std::to_string(c.flag);
            sam_out += "\t";
            if (c.flag & 4) { // unmapped here
                // RNAME POS MAPQ(0) * RNEXT PNEXT TLEN
                sam_out += "*\t0\t0\t*\t*\t0\t0\t";
            } else {
                std::string s2_name(hdr->target_name[c.tid]);
                std::string s1_name(l->get_other_name(s2_name));
                sam_out += s1_name.data(); // REF
                sam_out += "\t";
                sam_out += std::to_string(l->s2_to_s1(s2_name, c.pos, chroms_not_found) + 1); // POS
                sam_out += "\t";
                sam_out += std::to_string(c.qual); // QUAL
                sam_out += "\t";
                sam_out += l->cigar_s2_to_s1(s2_name, aln).data(); // CIGAR
                sam_out += "\t*\t0\t0\t"; // RNEXT PNEXT TLEN
            }
            // get query sequence
            std::string query_seq("");
            uint8_t* seq = bam_get_seq(aln);
            for (auto i = 0; i < c.l_qseq; ++i) {
                query_seq += seq_nt16_str[bam_seqi(seq, i)];
            }
            sam_out += query_seq.data();
            sam_out += "\t";
            // get quality
            std::string qual_seq("");
            uint8_t* qual = bam_get_qual(aln);
            if (qual[0] == 255) qual_seq = "*";
            else {
                for (auto i = 0; i < c.l_qseq; ++i) {
                    qual_seq += (char) (qual[i] + 33);
                }
            }
            sam_out += qual_seq.data();
            sam_out += "\t";
            // TODO reconcile any tags that also need to be added.
            sam_out += tag_to_string(aln).data();
            sam_out += "\n";
            // mutex.lock();
            fprintf(out_sam_fp, "%s", sam_out.c_str());
            // mutex.unlock();
            // it = std::next(it);
        }
    }
    for (int i = 0; i < chunk_size; i++){
        bam_destroy1(aln_vec[i]);
    }
    aln_vec.clear();
}

void lift_run(lift_opts args) {
    // if "-l" not specified, then create a liftover
    lift::LiftMap l = [&]{
        if (args.lift_fname != "") {
            std::ifstream in(args.lift_fname, std::ios::binary);
            return lift::LiftMap(in);
        } else if (args.vcf_fname != "") {
            return lift::LiftMap(lift_from_vcf(args.vcf_fname, args.sample, args.haplotype, args.name_map, args.length_map));
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
    for (auto i = 0; i < hdr->n_targets; i++){
        // if a contig is not found in vcf/lft, print it as it was before liftover
        if (std::find(contig_names.begin(), contig_names.end(), hdr->target_name[i]) == contig_names.end()){
            fprintf(out_sam_fp, "@SQ\tSN:%s\tLN:%u\n", hdr->target_name[i], hdr->target_len[i]);
        }
    }
    /* This is only supported in the latest dev release of htslib as of July 22. 
     * add back in at next htslib release
    kstring_t s = KS_INITIALIZE;
    sam_hdr_find_line_id(hdr, "PG", NULL, NULL, &s);
    fprintf(out_sam_fp, "%s\n", s.s);
    ks_free(&s);
    */
    char* prev_pg = get_PG(hdr);
    fprintf(out_sam_fp, "%s\n", prev_pg);
    free(prev_pg);
    fprintf(out_sam_fp, "@PG\tID:liftover\tPN:liftover\tCL:\"%s\"\n", args.cmd.data());
    
    // store chromosomes found in SAM but not in lft
    // use a vector to avoid printing out the same warning msg multiple times
    std::vector<std::string> chroms_not_found;
    const int num_threads = args.threads;
    const int chunk_size = args.chunk_size;
    std::vector<std::thread> threads;
    std::mutex mutex;
    for (int j = 0; j < num_threads; j++){
        threads.push_back(
            std::thread(read_and_lift, &mutex, sam_fp, out_sam_fp, hdr, chunk_size, &l, &chroms_not_found)
        );
    }
    for (int j = 0; j < num_threads; j++){
        if(threads[j].joinable())
            threads[j].join();
    }
    threads.clear();
    fclose(out_sam_fp);
}

lift::LiftMap lift_from_vcf(std::string fname, 
                            std::string sample, 
                            std::string haplotype, 
                            NameMap names, LengthMap lengths) {
    if (fname == "" && sample == "") {
        fprintf(stderr, "vcf file name and sample name are required!! \n");
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
        {"haplotype", required_argument, 0, 'g'},
        {"threads", required_argument, 0, 't'},
        {"chunk_size", required_argument, 0, 'T'},
        {"verbose", no_argument, &args.verbose, 1},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "v:s:p:l:a:g:n:k:t:T:", long_options, &long_index)) != -1) {
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
            case 'k':
                args.length_map = parse_length_map(optarg);
                break;
            case 't':
                args.threads = atoi(optarg);
                break;
            case 'T':
                args.chunk_size = atoi(optarg);
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
        fprintf(stderr, "invalid haplotype %s\n", args.haplotype.c_str());
        exit(1);
    }
    if (!strcmp(argv[optind], "lift")) {
        lift_run(args);
    } else if (!strcmp(argv[optind], "serialize")) {
        serialize_run(args);
    }
    return 0;
}
