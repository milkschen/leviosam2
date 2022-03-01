/*
 * aln.hpp
 *
 * Perform dynamic programming to align sequences
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#ifndef ALN_HPP
#define ALN_HPP

#include <htslib/sam.h>
#include <vector>
#include "chain.hpp"
#include "ksw2.h"

namespace Aln {

std::vector<uint32_t> align_extd2(
    const char *tseq, const char *qseq, int sc_mch, int sc_mis,
    const int q, const int e, const int8_t q2, const int8_t e2, const int w,
    const int zdrop, const int end_bonus, const int flag
);

};

#endif
