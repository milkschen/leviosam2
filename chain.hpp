#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <iostream>
#include <regex>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>

namespace chain {
class Interval {
    public:
    Interval() {
        source = "";
        target = "";
        offset = 0;
        source_start = 0;
        source_end = 0;
        same_strand = true;
    }

    Interval(
        std::string s, std::string t, int so, int se, int o, bool ss
    ) {
        source = s;
        target = t;
        offset = o;
        source_start = so;
        source_end = se;
        same_strand = ss;
    }

    void debug_print_interval() {
        std::string ss = (same_strand)? "+" : "-";
        std::cerr << source << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
    }

    std::string source;
    std::string target;
    int offset;
    int source_start, source_end;
    bool same_strand; // true: "+"; false: "-"
};

class ChainFile {

    public:
    // ChainFile();
    //  {
    //     init_rs();
    // }

    ChainFile(std::string fname);

    void init_bitvectors(std::string source, int source_length);
    void sort_interval_map();
    void sort_intervals(std::string contig);
    void debug_print_interval_map();
    void debug_print_intervals(std::string contig);

    int get_start_rank(std::string contig, int pos) {return start_rs1_map[contig](pos);}
    int get_end_rank(std::string contig, int pos) {return end_rs1_map[contig](pos);}

    private:
    void init_rs() {
        for (auto& it: this->start_bv_map) {
            std::string contig = it.first;
            std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_start = this->start_map.find(contig);
            if (find_start == this->start_map.end()) {
                this->start_map[contig] = sdsl::sd_vector<>(it.second);
                this->start_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
            }
            sdsl::util::init_support(this->start_rs1_map[contig], &this->start_map[contig]);
        }
        for (auto& it: this->end_bv_map) {
            std::string contig = it.first;
            std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_start = this->end_map.find(contig);
            if (find_start == this->end_map.end()) {
                this->end_map[contig] = sdsl::sd_vector<>(it.second);
                this->end_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
            }
            sdsl::util::init_support(this->end_rs1_map[contig], &this->end_map[contig]);
        }
    }

    void parse_chain_line(
        std::string line, std::string &source, std::string &target,
        int &source_offset, int &target_offset, bool &same_strand);

    std::unordered_map<std::string, std::vector<chain::Interval>> interval_map;
    // std::vector<chain::Interval> intervals;
    std::unordered_map<std::string, sdsl::bit_vector> start_bv_map;
    std::unordered_map<std::string, sdsl::bit_vector> end_bv_map;
    std::unordered_map<std::string, sdsl::sd_vector<>> start_map;
    std::unordered_map<std::string, sdsl::sd_vector<>> end_map;
    std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> start_rs1_map;
    std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> end_rs1_map;
    // sdsl::sd_vector<> start;
    // sdsl::sd_vector<> end;
    // sdsl::sd_vector<>::rank_1_type start_rs1;
    // sdsl::sd_vector<>::rank_1_type end_rs1;
    // sdsl::sd_vector<>::select_0_type ins_sls0;
    // sdsl::sd_vector<>::select_0_type del_sls0;
    // sdsl::sd_vector<>::select_0_type snp_sls0;
};

ChainFile* chain_open(std::string fname);

using IntervalMap = std::unordered_map<std::string, chain::Interval>;
};

#endif
