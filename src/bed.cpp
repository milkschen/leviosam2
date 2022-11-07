/*
 * bed.cpp
 *
 * Classes and routines to lift over genomic intervals in the BED format
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "bed.hpp"

#include <iostream>
#include <vector>

#include "leviosam_utils.hpp"

namespace BedUtils {

Bed::Bed(){};

Bed::Bed(const std::string &fn) { init(fn); }

/* Initiate a BED object
 */
void Bed::init(const std::string &fn) {
    bed_fn = fn;
    is_valid = true;
    std::ifstream bed_f(fn);
    int cnt = 0;
    if (bed_f.is_open()) {
        std::string line;
        while (std::getline(bed_f, line)) {
            cnt += 1;
            if (add_interval(line) == false) {
                std::cerr << "[E::Bed] Failed to update BED record " << line
                          << "\n";
            }
        }
        bed_f.close();
    }
    auto contig_cnt = index();
    std::cerr << "[I::Bed] Read " << cnt << " BED records (" << contig_cnt
              << " contigs) from " << fn << "\n";
}

/* Index the interval trees for each contig and returns the number of contigs
 */
int Bed::index() {
    int contig_cnt = 0;
    for (auto &itr : intervals) {
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

/* Interval intersect query
 * Args:
 *   isec_threshold:
 *    - <0: return true if overlap >= 1bp
 *    - 1>= `isec_threshold` > 0: return true if overlap fraction >= `isec_threshold`
 *    - `isec_threshold` > 1: return true if overlap size >= `isec_threshold` (bp)
*/
bool Bed::intersect(const std::string &contig, const size_t &pos1,
                    const size_t &pos2, const float &isec_threshold) {
    if (!is_valid) return false;
    if (pos2 <= pos1) return false;
    if (intervals.find(contig) == intervals.end()) return false;
    std::vector<size_t> isec;
    bool ovlp_status = intervals[contig].overlap(pos1, pos2, isec);
    if (!ovlp_status) return false;
    // <=0 means any overlap is good
    if (isec_threshold <= 0) return (isec.size() > 0);

    int max_overlap = 0;
    for (auto& _o: isec) {
        if (_o > max_overlap) max_overlap = _o;
    }
    if (isec_threshold > 0 && isec_threshold <= 1)
        return (max_overlap >= (pos2 - pos1) * isec_threshold);

    if (isec_threshold > 1) 
        return (max_overlap >= isec_threshold);

    return false;
    // return (isec.size() > 0);
}

// Point intersect query
bool Bed::intersect(const std::string &contig, const size_t &pos) {
    return Bed::intersect(contig, pos, pos + 1);
}

BedMap Bed::get_intervals() { return intervals; }

std::string Bed::get_fn() { return bed_fn; }

};  // namespace BedUtils
