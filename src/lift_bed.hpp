/*
 * lift_bed.cpp
 *
 * Lifting over genomic intervals in the BED format
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef LIFT_BED_HPP
#define LIFT_BED_HPP

#include <iostream>

#include "chain.hpp"

struct lift_bed_opts {
    std::string cmd = "";
    std::string bed_fname = "";
    std::string chainmap_fname = "";
    std::string outpre = "";
    int allowed_gaps = 500;
    int verbose = 0;
};

void print_lift_bed_help_msg();

int lift_bed_record(chain::ChainMap& cmap, const std::string& line,
                    const int& allowed_gaps);
void lift_bed(lift_bed_opts& args);
int lift_bed_run(int argc, char** argv);

#endif
