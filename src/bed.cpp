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

Bed::Bed(const std::string &fn) {
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
    std::cerr << "Read " << cnt << " BED records from " << contig_cnt << " contigs.\n";
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
    std::vector<std::string> fields = LevioSamUtils::split_str(line, "\t");
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


bool Bed::intersect(const std::string &contig, const size_t &pos) {
    if (intervals.find(contig) == intervals.end())
        return false;
    std::vector<size_t> isec;
    intervals[contig].overlap(pos, pos+1, isec);
    return (isec.size() > 0);
}

BedMap Bed::get_intervals() {
    return intervals;
}

}; // namespace

// int main(int argc, char** argv) {
//     auto b = BedUtils::Bed("../testdata/test.bed");
//     std::cerr << "chr1:1000 " << b.intersect("chr1", 1000) << "\n";
//     std::cerr << "chr1:1500 " << b.intersect("chr1", 1500) << "\n";
//     std::cerr << "chr10:100 " << b.intersect("chr10", 100) << "\n";
// 
//     // IITree<size_t, bool> tree;
//     // tree.add(12, 34, true); // add an interval
//     // tree.add(0, 23, true);
//     // tree.add(34, 56, true);
//     // tree.index(); // index
//     // std::vector<size_t> a;
//     // tree.overlap(22, 25, a); // retrieve overlaps
//     // for (size_t i = 0; i < a.size(); ++i)
//     //     printf("%d\t%d\t%d\n", tree.start(a[i]), tree.end(a[i]), tree.data(a[i]));
// }
