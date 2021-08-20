#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <iostream>
#include <unordered_map>
#include <regex>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <htslib/sam.h>

namespace chain {
/* Chain interval object
 * Each interval is a gapless alignment between the source and the dest references.
 */
class Interval {
    public:
        Interval();
        Interval(std::string t, int32_t so, int32_t se, int32_t o, bool ss);
        Interval(std::ifstream& in);
        void debug_print_interval();
        // Save to stream
        size_t serialize(std::ofstream& out) const;
        // Load from stream
        void load(std::istream& in);

        std::string target;
        int32_t offset;
        int32_t source_start, source_end;
        bool strand; // true: "+"; false: "-"
};

using BitVectorMap = std::unordered_map<std::string, sdsl::bit_vector>;
using SdVectorMap = std::unordered_map<std::string, sdsl::sd_vector<>>;
using SdRank1Map = std::unordered_map<std::string, sdsl::sd_vector<>::rank_1_type>;
using IntervalMap = std::unordered_map<std::string, std::vector<Interval>>;

class ChainMap {
    public:
        ChainMap(): verbose(0), allowed_cigar_changes(0) {}
        ~ChainMap() {}
        ChainMap(std::string fname, int verbose, int allowed_cigar_changes);
        ChainMap(std::ifstream& in, int verbose, int allowed_cigar_changes);
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
        bool update_interval_indexes(
            const std::string contig, const int32_t pos,
            int32_t &sidx, int32_t &eidx
        );

        std::string lift_contig(std::string contig, size_t pos);
        void push_cigar(
            std::vector<uint32_t> &cigar, uint32_t len, uint16_t op,
            const bool no_reduct);
        void lift_cigar(const std::string &contig, bam1_t* aln);
        void lift_cigar(
            const std::string &contig, bam1_t* aln,
            int start_sidx, int end_sidx, int num_clipped
        );
        std::vector<uint32_t> lift_cigar_core(
            const std::string &contig, bam1_t* aln, int num_clipped,
            const int start_sidx, const int end_sidx
        );

        size_t lift_pos(std::string contig, size_t pos);
        int32_t get_num_clipped(
            const int32_t pos, const bool leftmost,
            const std::string &contig, int32_t &sidx, int32_t &eidx
        );
        bool lift_segment(
            bam1_t* aln, bam_hdr_t* hdr,
            bool first_seg, std::string &dest_contig
        );
        void lift_aln(
            bam1_t* aln,
            bam_hdr_t* hdr,
            std::string &dest_contig
        );

        void parse_chain_line(
            std::string line, std::string &source, std::string &target,
            int32_t &source_len, int32_t &source_offset,
            int32_t &target_offset, bool &strand,
            BitVectorMap &start_bv_map, BitVectorMap &end_bv_map
        );

        bam_hdr_t* bam_hdr_from_chainmap(samFile* sam_fp);

        size_t serialize(std::ofstream& out);
        void load(std::ifstream& in);

    private:
        void init_rs();
        void lift_cigar_core_one_run(
            bam1_t* aln,
            std::vector<uint32_t> &new_cigar,
            std::queue<std::tuple<int32_t, int32_t>> &break_points,
            uint32_t cigar_op_len,
            unsigned int cigar_op,
            const uint32_t qlen,
            int &tmp_gap,
            int &query_offset
        );
        void update_flag_unmap(bam1_core_t* c, const bool first_seg);

        const int verbose;
        const int allowed_cigar_changes;
        IntervalMap interval_map;
        SdVectorMap start_map;
        SdVectorMap end_map;
        SdRank1Map start_rs1_map;
        SdRank1Map end_rs1_map;
        std::unordered_map<std::string, int32_t> length_map;

        // Debug functions
        void debug_print_interval_queries(
            const bool first_seg, const bool leftmost,
            const std::string contig, const int32_t pos,
            const int32_t sidx, const int32_t eidx
        );
};

};

#endif
