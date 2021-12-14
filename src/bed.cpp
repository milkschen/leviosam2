/*
 * bed.cpp
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#include <iostream>
#include <vector>
#include "bed.hpp"
#include "leviosam_utils.hpp"

namespace BedUtils {


Bed::Bed() {};


Bed::Bed(const std::string &fn) {
    init(fn);
}


/* Initiate a BED object
 */
void Bed::init(const std::string& fn) {
    bed_fn = fn;
    is_valid = true;
    std::ifstream bed_f(fn);
    int cnt = 0;
    if (bed_f.is_open()) {
        std::string line;
        while (std::getline(bed_f, line)) {
            cnt += 1;
            // std::cerr << line << "\n";
            if (add_interval(line) == false) {
                std::cerr << "[E::Bed] Failed to update BED record " << line << "\n";
            }
        }
        bed_f.close();
    }
    auto contig_cnt = index();
    std::cerr << "[I::Bed] Read " << cnt << " BED records (" << contig_cnt << " contigs) from " << fn << "\n";
}


int Bed::index() {
    int contig_cnt = 0;
    // Index the interval trees for each contig
    for (auto& itr: intervals) {
        itr.second.index();
        contig_cnt += 1;
    }
    return contig_cnt;
}


bool Bed::add_interval(const std::string &line) {
    is_valid = true;
    std::vector<std::string> fields = LevioSamUtils::str_to_vector(line, "\t");
    if (fields.size() < 3) {
        return false;
    }
    std::string contig = fields[0];
    size_t start = std::stoi(fields[1]);
    size_t end = std::stoi(fields[2]);
    auto itr = intervals.find(contig);
    if (itr == intervals.end()) {
        IITree<size_t, bool> new_tree;
        new_tree.add(start, end, true);
        intervals.emplace(std::make_pair(contig, new_tree));
    } else {
        itr->second.add(start, end, true);
    }
    return true;
}


// Interval intersect query
bool Bed::intersect(
    const std::string &contig, const size_t &pos1, const size_t &pos2) {
    if (!is_valid)
        return false;
    if (intervals.find(contig) == intervals.end())
        return false;
    std::vector<size_t> isec;
    intervals[contig].overlap(pos1, pos2, isec);
    return (isec.size() > 0);
}


// Point intersect query
bool Bed::intersect(const std::string &contig, const size_t &pos) {
    return Bed::intersect(contig, pos, pos+1);
}


BedMap Bed::get_intervals() {
    return intervals;
}


std::string Bed::get_fn() {
    return bed_fn;
}


}; // namespace

