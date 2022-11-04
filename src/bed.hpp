/*
 * bed.hpp
 *
 * Classes and routines to lift over genomic intervals in the BED format
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef BED_HPP
#define BED_HPP

#include <iostream>
#include <unordered_map>

#include "IITree.h"
#include "robin_hood.h"

namespace BedUtils {

using BedMap =
    robin_hood::unordered_map<std::string, IITree<std::size_t, bool>>;

class Bed {
   public:
    Bed();
    Bed(const std::string &fn);
    void init(const std::string &fn);

    int index();
    bool add_interval(const std::string &line);
    bool intersect(const std::string &contig, const size_t &pos1,
                   const size_t &pos2, const float &isec_frac = 0);
    bool intersect(const std::string &contig, const size_t &pos);
    BedMap get_intervals();
    std::string get_fn();

   private:
    BedMap intervals;
    std::string bed_fn = "";
    bool is_valid = false;

};  // Bed class

};  // namespace BedUtils

#endif
