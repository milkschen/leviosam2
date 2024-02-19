/*
 * cigar.hpp
 *
 * Cigar-related utilities.
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef CIGAR_HPP
#define CIGAR_HPP

#include <htslib/sam.h>

#include <iostream>
#include <vector>

#include "leviosam_utils.hpp"

namespace Cigar {

using CigarVector = std::vector<uint32_t>;

/// Adds a cigar element to a cigar vector.
void push_cigar(CigarVector& cigar_vec, uint32_t len, uint16_t op,
                const bool no_reduce);
void pop_cigar(CigarVector& cigar_vec, uint32_t size);
void sclip_cigar_front(uint32_t* cigar, const uint32_t& n_cigar, int len_clip,
                       CigarVector& new_cigar_vec, int& idx, int& query_offset);
void sclip_cigar_back(CigarVector& cigar_vec, int len_clip);
void set_empty_cigar(bam1_t* aln);
void update_cigar(bam1_t* aln, CigarVector& new_cigar_vec);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);

}  // namespace Cigar

#endif