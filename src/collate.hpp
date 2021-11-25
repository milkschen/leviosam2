/*
 * collate.hpp
 *
 * The `leviosam collate` program that makes sure a split paired-end
 * BAM file is properly paired
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#ifndef COLLATE_HPP
#define COLLATE_HPP
#include <string>


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

void print_collate_help_msg();

int collate_run(int argc, char** argv);

#endif
