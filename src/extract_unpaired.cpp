
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include <htslib/sam.h>
#include "leviosam.hpp"
#include "leviosam_utils.cpp"
#include "robin_hood.h"


struct extract_unpaired_opts {
    std::string cmd = "";
    std::string sam_fname = "";
    std::string deferred_sam_fname = "";
    std::string outpre = "";
    std::string fq_fname = "";

    // Non-arguments
    std::string out_deferred_sam_fname = "";
    std::string out_committed_sam_fname = "";
    std::string out_r1_fname = "";
    std::string out_r2_fname = "";
};

std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}


void print_help_msg() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: extract_unpaired (extract unpaired alignments from a BAM using a FASTQ of singletons)\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Usage:   extract_unpaired [options] -a <bam> {-b <bam> | -q <fastq>} -o <prefix>\n\n");
    fprintf(stderr, "Inputs:  -a string   Path to the input SAM/BAM.\n");
    fprintf(stderr, "         -b string   Path to the input deferred SAM/BAM.\n");
    fprintf(stderr, "         -q string   Path to the input singleton FASTQ.\n");
    fprintf(stderr, "         -o string   Prefix to the output files (1 BAM, a pair of FASTQs).\n");
    fprintf(stderr, "Options: -h          Print detailed usage.\n");
    fprintf(stderr, "         -V INT      Verbose level [0].\n");
    fprintf(stderr, "\n");
}


void extract_unpaired_core(
    LevioSamUtils::fastq_map &reads,
    bam_hdr_t* chdr,
    bam_hdr_t* dhdr,
    samFile* csam_fp, samFile* out_csam_fp,
    samFile* out_dsam_fp,
    std::ofstream& out_r1_fp, std::ofstream& out_r2_fp
) {
    bam1_t* aln = bam_init1();
    int cnt = 0;
    while (sam_read1(csam_fp, chdr, aln) > 0) {
        bam1_core_t c = aln->core;
        if ((c.flag & 256) || // Secondary alignment - no SEQ field
            (c.flag & 512) || // not passing filters
            (c.flag & 1024) || // PCR or optinal duplicate
            (c.flag & 2048)) { // supplementary alignment
            if (sam_write1(out_csam_fp, chdr, aln) < 0) {
                std::cerr << "[Error] Failed to write record " << 
                    bam_get_qname(aln) << "\n";
                exit(1);
            }
        } 
        std::string qname = bam_get_qname(aln);
        auto search = reads.find(qname);
        bool write_to_fastq = false;
        if (search != reads.end()){
            LevioSamUtils::FastqRecord fq = LevioSamUtils::FastqRecord(aln);
            // Write the paired records to FASTQ
            if (c.flag & 64) { // first segment
                fq.write(out_r1_fp, qname);
                search->second.write(out_r2_fp, qname);
            } else if (c.flag & 128) { // second segment
                search->second.write(out_r1_fp, qname);
                fq.write(out_r2_fp, qname);
            } else {
                std::cerr << "Error: Read " << qname << " is not paired-end. Exit.\n";
                exit(1);
            }
            // Also write the records to a BAM file
            if (sam_write1(out_dsam_fp, dhdr, aln) < 0 ||
                sam_write1(out_dsam_fp, dhdr, search->second.aln) < 0) {
                std::cerr << "[Error] Failed to write record " << bam_get_qname(aln) <<
                    " to the deferred BAM file\n";
                exit(1);
            }
            reads.erase(search);
            write_to_fastq = true;
            cnt += 1;
        }
        if (write_to_fastq == false) { // write to BAM
            if (sam_write1(out_csam_fp, chdr, aln) < 0) {
                std::cerr << "[Error] Failed to write record " << bam_get_qname(aln) << 
                    " to the committed BAM file\n";
                exit(1);
            }
        }
    }
    
    std::cerr << "Extract " << cnt << " reads from BAM\n";
    std::cerr << "Num. remaining records in the map = " << reads.size() << " (expected to be 0)\n";
}


void extract_unpaired(extract_unpaired_opts args) {
    // Input file
    samFile* csam_fp = (args.sam_fname == "")?
        sam_open("-", "r") : sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* chdr = sam_hdr_read(csam_fp);
    samFile* dsam_fp = (args.deferred_sam_fname == "")?
        NULL : sam_open(args.deferred_sam_fname.data(), "r");
    bam_hdr_t* dhdr = (args.deferred_sam_fname == "")?
        NULL : sam_hdr_read(dsam_fp);
    // Output files
    std::ofstream out_r1_fp(args.out_r1_fname);
    std::ofstream out_r2_fp(args.out_r2_fname);
    samFile* out_csam_fp = sam_open(args.out_committed_sam_fname.data(), "wb");
    samFile* out_dsam_fp = sam_open(args.out_deferred_sam_fname.data(), "wb");

    sam_hdr_add_pg(chdr, "extract_unpaired", "VN", VERSION, "CL", args.cmd.data(), NULL);
    sam_hdr_add_pg(dhdr, "extract_unpaired", "VN", VERSION, "CL", args.cmd.data(), NULL);
    if (sam_hdr_write(out_csam_fp, chdr) < 0 || sam_hdr_write(out_dsam_fp, dhdr) < 0) {
        std::cerr << "Error: Unable to write SAM header\n";
        exit(1);
    }

    // Core operation
    LevioSamUtils::fastq_map reads = (args.fq_fname != "")?
        LevioSamUtils::read_unpaired_fq(args.fq_fname) :
        LevioSamUtils::read_deferred_bam(dsam_fp, out_dsam_fp, dhdr, out_r1_fp, out_r2_fp);
    extract_unpaired_core(
        reads, chdr, dhdr, csam_fp, out_csam_fp, out_dsam_fp, out_r1_fp, out_r2_fp);

    if (dsam_fp != NULL)
        sam_close(dsam_fp);
    sam_close(out_dsam_fp);
    out_r1_fp.close();
    out_r2_fp.close();
    sam_close(csam_fp);
    sam_close(out_csam_fp);
}


int main(int argc, char** argv) {
    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();
    int c;
    extract_unpaired_opts args;
    args.cmd = make_cmd(argc,argv);
    static struct option long_options[] {
        {"sam", required_argument, 0, 'a'},
        {"deferred_sam", required_argument, 0, 'b'},
        {"output", required_argument, 0, 'p'},
        {"fastq", required_argument, 0, 'q'},
        {"verbose", required_argument, 0, 'V'},
    };
    int long_index = 0;
    while((c = getopt_long(argc, argv, "ha:b:p:q:V:", long_options, &long_index)) != -1) {
        switch (c) {
            case 'h':
                print_help_msg();
                exit(0);
            case 'a':
                args.sam_fname = optarg;
                break;
            case 'b':
                args.deferred_sam_fname = optarg;
                break;
            case 'p':
                args.outpre = optarg;
                break;
            case 'q':
                args.fq_fname = optarg;
                break;
            default:
                fprintf(stderr, "ignoring option %c\n", c);
                exit(1);
        }
    }
    if (args.sam_fname == "" || args.outpre == "" ||
        (args.fq_fname == "" && args.deferred_sam_fname == "")) {
        std::cerr << "Error: required argument missed.\n";
        print_help_msg();
        exit(1);
    }
    if (args.fq_fname != "" && args.deferred_sam_fname != "") {
        std::cerr << "Error: only one of `-q` and `-b` can be set.\n";
        print_help_msg();
        exit(1);
    }

    std::cerr << "Inputs:\n";
    std::string input_sam = (args.sam_fname == "" || args.sam_fname == "-")? "stdin" : args.sam_fname;
    std::cerr << " - BAM: " << input_sam << "\n";
    if (args.deferred_sam_fname != "")
        std::cerr << " - BAM (deferred): " << args.deferred_sam_fname << "\n";
    if (args.fq_fname != "")
        std::cerr << " - FASTQ: " << args.fq_fname << "\n";

    std::cerr << "\nOutputs:\n";
    args.out_committed_sam_fname = args.outpre + "-committed.bam";
    std::cerr << " - BAM (committed): " << args.out_committed_sam_fname << "\n";
    if (args.deferred_sam_fname != "") {
        args.out_deferred_sam_fname = args.outpre + "-deferred.bam";
        std::cerr << " - BAM (deferred): " << args.out_deferred_sam_fname << "\n";
    }

    args.out_r1_fname = args.outpre + "-deferred-R1.fq";
    args.out_r2_fname = args.outpre + "-deferred-R2.fq";
    std::cerr << " - FASTQ1: " << args.out_r1_fname << "\n";
    std::cerr << " - FASTQ2: " << args.out_r2_fname + "\n";
    std::cerr << "\n";

    extract_unpaired(args);
    
    double cpu_duration = (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration = (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "Finished in " << cpu_duration << " CPU seconds, or " << 
                                   wall_duration.count() << " wall clock seconds\n";

    return 0;
}
