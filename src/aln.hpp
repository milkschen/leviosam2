/*
 * aln.hpp
 *
 * Perform dynamic programming to align sequences
 *
 * Author: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef ALN_HPP
#define ALN_HPP

#include <htslib/sam.h>

#include <vector>

#include "chain.hpp"
#include "cigar.hpp"
#include "ksw2.h"
#include "leviosam_utils.hpp"
#include "yaml.hpp"

namespace Aln {

class AlnOpts {
   public:
    std::string engine = "ksw_extz";
    int flag = 0;
    int nm_threshold = 15;
    int a = 2, b = 4, q = 4, e = 2, q2 = 24, e2 = 1;
    int w = 601, zdrop = 400, end_bonus = -1;

    void deserialize_realn(ryml::Tree realign_tree);
    void print_parameters();
};

int align_ksw2(const char* tseq, const char* qseq, const AlnOpts& opt,
               std::vector<uint32_t>& new_cigar, int& new_score);

};  // namespace Aln

#endif
