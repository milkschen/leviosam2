/*
 * collate.cpp
 *
 * The `leviosam collate` program that collates a pair of BAM files that are originally split
 * from a paired-end BAM file (example use case: a BAM is split using a MAPQ cutoff). The
 * resulting pair of files will be properly paired
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#include <ctime>
#include <getopt.h>
#include <iostream>
#include <stdio.h>
#include "collate.hpp"
#include "version.hpp"


void print_collate_help_msg() {
    fprintf(stderr, "\n");
    fprintf(stderr, "Collate alignments to make sure reads are paired\n");
    fprintf(stderr, "Version: %s\n", VERSION);
    fprintf(stderr, "Usage:   leviosam collate [options] -a <bam> {-b <bam> | -q <fastq>} -p <prefix>\n\n");
    fprintf(stderr, "Inputs:  -a string   Path to the input SAM/BAM.\n");
    fprintf(stderr, "         -b string   Path to the input deferred SAM/BAM.\n");
    fprintf(stderr, "         -q string   Path to the input singleton FASTQ.\n");
    fprintf(stderr, "         -p string   Prefix to the output files (1 BAM and a pair of gzipped FASTQs).\n");
    fprintf(stderr, "Options: -h          Print detailed usage.\n");
    fprintf(stderr, "         -V INT      Verbose level [0].\n");
    fprintf(stderr, "\n");
}


/* Read a SAM/BAM file, write properly paired alignements to a BAM file
 * and return the rest as a fastq_map
 */
fastq_map read_deferred_bam(
    samFile* dsam_fp, samFile* out_dsam_fp, sam_hdr_t* hdr,
    ogzstream& out_r1_fp, ogzstream& out_r2_fp
) {
    fastq_map reads1, reads2;
    bam1_t* aln = bam_init1();

    while (sam_read1(dsam_fp, hdr, aln) > 0) {
        bam1_core_t c = aln->core;
        // The following categories of reads are excluded by this method
        if ((c.flag & BAM_FSECONDARY) || // Secondary alignment - no SEQ field
            (c.flag & BAM_FQCFAIL) || // not passing filters
            (c.flag & BAM_FDUP) || // PCR or optinal duplicate
            (c.flag & BAM_FSUPPLEMENTARY)) { // supplementary alignment
            continue;
        }
        std::string qname = bam_get_qname(aln);
        LevioSamUtils::FastqRecord fq = LevioSamUtils::FastqRecord(aln);
        if (c.flag & BAM_FREAD1) { // first segment
            auto search = reads2.find(qname);
            if (search != reads2.end()) {
                if (sam_write1(out_dsam_fp, hdr, search->second.aln) < 0  ||
                    sam_write1(out_dsam_fp, hdr, aln) < 0) {
                    std::cerr << "[Error] Failed to write record " << qname << "\n";
                    exit(1);
                }
                fq.write(out_r1_fp, qname);
                search->second.write(out_r2_fp, qname);
                reads2.erase(search);
            } else {
                reads1.emplace(std::make_pair(qname, fq));
            }
        } else if (c.flag & BAM_FREAD2) { // second segment
            auto search = reads1.find(qname);
            if (search != reads1.end()) {
                if (sam_write1(out_dsam_fp, hdr, aln) < 0 ||
                    sam_write1(out_dsam_fp, hdr, search->second.aln) < 0) {
                    std::cerr << "[Error] Failed to write record " << qname << "\n";
                    exit(1);
                }
                search->second.write(out_r1_fp, qname);
                fq.write(out_r2_fp, qname);
                reads1.erase(search);
            } else {
                reads2.emplace(std::make_pair(qname, fq));
            }
        } else {
            std::cerr << "[Error] Alignment " << qname << " is not paired.\n";
            exit(1);
        }
    }
    std::cerr << "Num. unpaired reads1: " << reads1.size() << "\n";
    std::cerr << "Num. unpaired reads2: " << reads2.size() << "\n";
    for (auto& r: reads2) {
        reads1.emplace(std::make_pair(r.first, r.second.aln));
    }
    std::cerr << "Num. merged unpaired: " << reads1.size() << "\n";
    return reads1;
}


/* Read a single-end FASTQ file and return a fastq_map */
fastq_map read_unpaired_fq(const std::string& fq_fname) {
    fastq_map reads;
    std::ifstream fastq_fp(fq_fname);
    std::string line;
    int i = 0;
    std::string name, seq;
    while (getline (fastq_fp, line)) {
        if (i % 4 == 0) {
            name = line.substr(1);
        } else if (i % 4 == 1) {
            seq = line;
        } else if (i % 4 == 3) {
            reads.emplace(
                std::make_pair(name, LevioSamUtils::FastqRecord(seq, line)));
            name = "";
            seq = "";
        }
        i++;
    }
    std::cerr << "Number of singletons: " << reads.size() << "\n";
    fastq_fp.close();
    return reads;
}


void collate_core(
    fastq_map &reads,
    bam_hdr_t* chdr,
    bam_hdr_t* dhdr,
    samFile* csam_fp, samFile* out_csam_fp,
    samFile* out_dsam_fp,
    ogzstream& out_r1_fp, ogzstream& out_r2_fp
) {
    bam1_t* aln = bam_init1();
    int cnt = 0;
    while (sam_read1(csam_fp, chdr, aln) > 0) {
        bam1_core_t c = aln->core;
        bool primary = true;
        if ((c.flag & BAM_FSECONDARY) || // Secondary alignment - no SEQ field
            (c.flag & BAM_FQCFAIL) || // not passing filters
            (c.flag & BAM_FDUP) || // PCR or optinal duplicate
            (c.flag & BAM_FSUPPLEMENTARY)) { // supplementary alignment
            primary = false;
            // if (sam_write1(out_csam_fp, chdr, aln) < 0) {
            //     std::cerr << "[Error] Failed to write record " << 
            //         bam_get_qname(aln) << "\n";
            //     exit(1);
            // }
        } 
        std::string qname = bam_get_qname(aln);
        auto search = reads.find(qname);
        bool write_to_fastq = false;
        if (search != reads.end()){
            // If a read is not primary but in the defer group, skip it
            if (!primary)
                continue;

            LevioSamUtils::FastqRecord fq = LevioSamUtils::FastqRecord(aln);
            // Write the paired records to FASTQ
            if (c.flag & BAM_FREAD1) { // first segment
                fq.write(out_r1_fp, qname);
                search->second.write(out_r2_fp, qname);
            } else if (c.flag & BAM_FREAD2) { // second segment
                search->second.write(out_r1_fp, qname);
                fq.write(out_r2_fp, qname);
            } else {
                std::cerr << "[E::collate_core] Read " << qname << " is not paired-end. Exit.\n";
                exit(1);
            }
            // Also write the records to a BAM file
            if (sam_write1(out_dsam_fp, dhdr, aln) < 0 ||
                sam_write1(out_dsam_fp, dhdr, search->second.aln) < 0) {
                std::cerr << "[E::collate_core] Failed to write record " << bam_get_qname(aln) <<
                    " to the deferred BAM file\n";
                exit(1);
            }
            reads.erase(search);
            write_to_fastq = true;
            cnt += 1;
        } else { // read not found in defer, write to committed BAM
        // if (write_to_fastq == false) { // write to BAM
            if (sam_write1(out_csam_fp, chdr, aln) < 0) {
                std::cerr << "[E::collate_core] Failed to write record " << bam_get_qname(aln) << 
                    " to the committed BAM file\n";
                exit(1);
            }
        }
    }
    
    std::cerr << "[I::collate_core] Extract " << cnt << " reads from BAM\n";
    if (reads.size() != 0) {
        std::cerr << "[W::collate_core] Num. remaining records in the map = "
                  << reads.size() << " (expected to be 0)\n";
    }
}


void collate(collate_opts args) {
    // Input file
    samFile* csam_fp = (args.sam_fname == "")?
        sam_open("-", "r") : sam_open(args.sam_fname.data(), "r");
    bam_hdr_t* chdr = sam_hdr_read(csam_fp);
    samFile* dsam_fp = (args.deferred_sam_fname == "")?
        NULL : sam_open(args.deferred_sam_fname.data(), "r");
    bam_hdr_t* dhdr = (args.deferred_sam_fname == "")?
        NULL : sam_hdr_read(dsam_fp);

    // Output files
    ogzstream out_r1_fp(args.out_r1_fname.data());
    ogzstream out_r2_fp(args.out_r2_fname.data());
    samFile* out_csam_fp = sam_open(args.out_committed_sam_fname.data(), "wb");
    samFile* out_dsam_fp = sam_open(args.out_deferred_sam_fname.data(), "wb");

    sam_hdr_add_pg(chdr, "leviosam", "VN", VERSION, "CL", args.cmd.data(), NULL);
    sam_hdr_add_pg(dhdr, "leviosam", "VN", VERSION, "CL", args.cmd.data(), NULL);
    if (sam_hdr_write(out_csam_fp, chdr) < 0 || sam_hdr_write(out_dsam_fp, dhdr) < 0) {
        std::cerr << "Error: Unable to write SAM header\n";
        exit(1);
    }

    // Core operation
    fastq_map reads = (args.fq_fname != "")?
        read_unpaired_fq(args.fq_fname) :
        read_deferred_bam(dsam_fp, out_dsam_fp, dhdr, out_r1_fp, out_r2_fp);
    collate_core(
        reads, chdr, dhdr, csam_fp, out_csam_fp, out_dsam_fp, out_r1_fp, out_r2_fp);

    if (dsam_fp != NULL)
        sam_close(dsam_fp);
    sam_close(out_dsam_fp);
    out_r1_fp.close();
    out_r2_fp.close();
    sam_close(csam_fp);
    sam_close(out_csam_fp);
}


int collate_run(int argc, char** argv) {
    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();
    int c;
    collate_opts args;
    args.cmd = LevioSamUtils::make_cmd(argc,argv);
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
                print_collate_help_msg();
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
        std::cerr << "[E::collate_run] required argument missed.\n";
        print_collate_help_msg();
        exit(1);
    }
    if (args.fq_fname != "" && args.deferred_sam_fname != "") {
        std::cerr << "[E::collate_run] only one of `-q` and `-b` can be set.\n";
        print_collate_help_msg();
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

    args.out_r1_fname = args.outpre + "-deferred-R1.fq.gz";
    args.out_r2_fname = args.outpre + "-deferred-R2.fq.gz";
    std::cerr << " - FASTQ1: " << args.out_r1_fname << "\n";
    std::cerr << " - FASTQ2: " << args.out_r2_fname + "\n";
    std::cerr << "\n";

    collate(args);
    
    double cpu_duration = (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration = (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "Finished in " << cpu_duration << " CPU seconds, or " << 
                                   wall_duration.count() << " wall clock seconds\n";

    return 0;
}
