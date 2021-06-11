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
        // source = "";
        target = "";
        offset = 0;
        source_start = 0;
        source_end = 0;
        same_strand = true;
    }

    Interval(
        // std::string s, 
        std::string t, int so, int se, int o, bool ss
    ) {
        // source = s;
        target = t;
        offset = o;
        source_start = so;
        source_end = se;
        same_strand = ss;
    }

    void debug_print_interval() {
        std::string ss = (same_strand)? "+" : "-";
        std::cerr << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
        // std::cerr << source << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
    }

    // std::string source;
    std::string target;
    int offset;
    int source_start, source_end;
    bool same_strand; // true: "+"; false: "-"
};

class ChainFile {

    public:
    ChainFile(std::string fname, int verbose);
    void init_bitvectors(std::string source, int source_length);
    void sort_interval_map();
    void sort_intervals(std::string contig);
    void debug_print_interval_map();
    void debug_print_intervals(std::string contig);
    bool interval_map_sanity_check();
    int get_start_rank(std::string contig, int pos);
    int get_end_rank(std::string contig, int pos);
    void show_interval_info(std::string contig, int pos);
    int get_lifted_pos(std::string contig, int pos);

    private:
    void init_rs();

    void parse_chain_line(
        std::string line, std::string &source, std::string &target,
        int &source_offset, int &target_offset, bool &same_strand);

    int verbose;
    std::unordered_map<std::string, std::vector<chain::Interval>> interval_map;
    std::unordered_map<std::string, sdsl::bit_vector> start_bv_map;
    std::unordered_map<std::string, sdsl::bit_vector> end_bv_map;
    std::unordered_map<std::string, sdsl::sd_vector<>> start_map;
    std::unordered_map<std::string, sdsl::sd_vector<>> end_map;
    std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> start_rs1_map;
    std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> end_rs1_map;
};

ChainFile* chain_open(std::string fname, int verbose);

using IntervalMap = std::unordered_map<std::string, chain::Interval>;
};

#endif
