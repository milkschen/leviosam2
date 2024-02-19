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

namespace Cigar {

void push_cigar(std::vector<uint32_t>& cigar, uint32_t len, uint16_t op,
                const bool no_reduct);
void pop_cigar(std::vector<uint32_t>& cigar, uint32_t size);
void sclip_cigar_front(uint32_t* cigar, const uint32_t& n_cigar, int len_clip,
                       std::vector<uint32_t>& new_cigar, int& idx,
                       int& query_offset);
void sclip_cigar_back(std::vector<uint32_t>& new_cigar, int len_clip);

void update_cigar(bam1_t* aln, std::vector<uint32_t>& new_cigar);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);

}  // namespace Cigar

#endif