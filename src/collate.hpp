/*
 * collate.hpp
 *
 * Collates a pair of BAM files that are originally split from a paired-end
 * dataset. The resulting pair of files will be properly paired
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef COLLATE_HPP
#define COLLATE_HPP
#include <htslib/sam.h>

#include <string>

#include "gzstream.h"
#include "leviosam_utils.hpp"
#include "robin_hood.h"

struct collate_opts {
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

typedef robin_hood::unordered_map<std::string, LevioSamUtils::FastqRecord>
    fastq_map;

fastq_map read_unpaired_fq(const std::string& fq_fname);
fastq_map read_deferred_bam(samFile* dsam_fp, samFile* out_dsam_fp,
                            sam_hdr_t* hdr, ogzstream& out_r1_fp,
                            ogzstream& out_r2_fp);

void print_collate_help_msg();

int collate_run(int argc, char** argv);

#endif
