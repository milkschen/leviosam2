#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <iostream>
#include <regex>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <htslib/sam.h>
#include "leviosam.hpp"


namespace chain {
class Interval {
    public:
    Interval() {
        // source = "";
        target = "";
        offset = 0;
        source_start = 0;
        source_end = 0;
        strand = true;
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
        strand = ss;
    }
    
    Interval(std::ifstream& in) {
        this->load(in);
    }

    void debug_print_interval() {
        std::string ss = (strand)? "+" : "-";
        std::cerr << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
        // std::cerr << source << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
    }

    // std::string source;
    std::string target;
    int offset;
    int source_start, source_end;
    bool strand; // true: "+"; false: "-"


    // Save to stream
    size_t serialize(std::ofstream& out) const;
    // Load from stream
    void load(std::istream& in);

    void debug_print();
};

using BitVectorMap = std::unordered_map<std::string, sdsl::bit_vector>;
using SdVectorMap = std::unordered_map<std::string, sdsl::sd_vector<>>;
using IntervalMap = std::unordered_map<std::string, std::vector<Interval>>;


class ChainMap {

    public:
        ChainMap() {}
        ~ChainMap() {}
        ChainMap(std::string fname, int verbose);
        ChainMap(std::ifstream& in);
        void init_bitvectors(
            std::string source, int source_length,
            std::unordered_map<std::string, sdsl::bit_vector> &start_bv_map,
            std::unordered_map<std::string, sdsl::bit_vector> &end_bv_map);
        void sort_interval_map();
        void sort_intervals(std::string contig);
        void debug_print_interval_map();
        void debug_print_intervals(std::string contig);
        bool interval_map_sanity_check();
        int get_start_rank(std::string contig, int pos);
        int get_end_rank(std::string contig, int pos);

        void show_interval_info(std::string contig, int pos);

        std::string lift_contig(std::string contig, size_t pos);
        std::string lift_contig(
            std::string contig, int start_intvl_idx, int end_intvl_idx);
        void lift_cigar(const std::string& contig, bam1_t* aln);
        void lift_cigar(
            const std::string& contig, bam1_t* aln,
            int start_intvl_idx, int end_intvl_idx);

        size_t lift_pos(std::string contig, size_t pos);
        size_t lift_pos(
            std::string contig, size_t pos,
            int start_intvl_idx, int end_intvl_idx);

        bool lift_segment(
            bam1_t* aln, bam_hdr_t* hdr,
            bool is_first_seg, std::string &dest_contig,
            int &start_intvl_idx, int &end_intvl_idx);
        bool is_liftable(
            std::string contig, size_t pos,
            int &start_intvl_idx, int &end_intvl_idx);
        void lift_aln(
            bam1_t* aln,
            bam_hdr_t* hdr,
            std::string &dest_contig);
            // bool md_flag,
            // std::string &ref_name,
            // std::map<std::string, std::string>* ref_dict);


        void parse_chain_line(
            std::string line, std::string &source, std::string &target,
            int &source_len, int &source_offset, int &target_offset, bool &strand,
            BitVectorMap &start_bv_map, BitVectorMap &end_bv_map);

        size_t serialize(std::ofstream& out);
        void load(std::ifstream& in);

    private:
        void init_rs();

        int verbose;
        IntervalMap interval_map;
        SdVectorMap start_map;
        SdVectorMap end_map;
        std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> start_rs1_map;
        std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type> end_rs1_map;
};

};

#endif
