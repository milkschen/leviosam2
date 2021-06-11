#include <iostream>
#include "chain.hpp"

using namespace chain;

ChainFile::ChainFile(std::string fname, int verbose) {
    this->verbose = verbose;
    std::ifstream chain_f(fname);
    std::string line;
    if (chain_f.is_open())
    {
        std::string source;
        std::string target;
        int source_offset = 0;
        int target_offset = 0;
        bool current_ss = true;
        while (getline(chain_f, line)) {
            this->parse_chain_line(
                line, source, target, source_offset, 
                target_offset, current_ss);
        }
        chain_f.close();
    }
    this->init_rs();
    if (this->verbose > 1) {
        this->debug_print_interval_map();
        std::cerr << "\n\nsort intervals\n\n";
    }
    this->sort_interval_map();
    if (this->verbose > 1) {
        this->debug_print_interval_map();
    }
}


/* Create start and end bitvectors when see a new `source`.
 * Do nothing if `source` has been seen.
 */
void ChainFile::init_bitvectors(std::string source, int source_length) {
    std::unordered_map<std::string, sdsl::bit_vector>::const_iterator find_start = this->start_bv_map.find(source);
    if (find_start == this->start_bv_map.end()) {
        sdsl::bit_vector new_start_bv(source_length), new_end_bv(source_length);
        this->start_bv_map[source] = new_start_bv;
        this->end_bv_map[source] = new_end_bv;
    }
}


void ChainFile::sort_interval_map() {
    for (auto &itr: this->interval_map) {
        this->sort_intervals(itr.first);
    }
}


/* Sort the intervals in a ChainFile in ascending order. */
void ChainFile::sort_intervals(std::string contig) {
    std::sort(
        this->interval_map[contig].begin(), this->interval_map[contig].end(),
        [](const Interval& lhs, const Interval& rhs) {
            return lhs.source_start < rhs.source_start;}
    );
}


void ChainFile::debug_print_interval_map() {
    for (auto &intvl_map: this->interval_map) {
        this->debug_print_intervals(intvl_map.first);
    }
}


/* Print the intervals in a ChainFile. */
void ChainFile::debug_print_intervals(std::string contig) {
    for (auto &intvl: this->interval_map[contig]) {
        intvl.debug_print_interval();
    }
    std::cerr << "\n";
}


/* Check if the interval map contains any overlaps.
 * Logic:
 *   For each interval, its ending position <= next starting position
 * 
 * Returns a bool:
 *   true if pass; false otherwise
 */
bool ChainFile::interval_map_sanity_check() {
    for (auto& itr : this->interval_map){
            std::vector<chain::Interval> v = this->interval_map[itr.first];
            for (int i = 0; i < v.size()-1; i++) {
                if (v[i].source_end > v[i+1].source_start) {
                    std::cerr << "Error: " << itr.first << "\n";
                    v[i].debug_print_interval();
                    return false;
                }
            }
    }
    return true;
}

/* Get the rank in the start bitvector at contig[pos]
 * If pos > len(contig): return rank(contig[-1])
 * Returns:
 *   non-negative int: rank
 *   -1: if pos < 0 or contig does not exist
 */
int ChainFile::get_start_rank(std::string contig, int pos) {
    if (pos < 0) return -1;
    std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_start = 
        this->start_map.find(contig);
    if (find_start == this->start_map.end()) {
        return -1;
    }
    if (pos >= start_rs1_map[contig].size())
        pos = start_rs1_map[contig].size() - 1;
    return start_rs1_map[contig](pos) - 1;
}

/* Get the rank in the end bitvector at contig[pos]
 * If pos > len(contig): return rank(contig[-1])
 * Returns:
 *   non-negative int: rank
 *   -1: if pos < 0 or contig does not exist
 */
int ChainFile::get_end_rank(std::string contig, int pos) {
    if (pos < 0) return -1;
    std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_end = this->end_map.find(contig);
    if (find_end == this->end_map.end()) {
        return -1;
    }
    if (pos >= end_rs1_map[contig].size())
        pos = end_rs1_map[contig].size() - 1;
    return end_rs1_map[contig](pos) - 1;
}

/* Init rank support for all start/end bitvectors. */
void ChainFile::init_rs() {
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

/* Parse a chain line and update the chain object.
 */
void ChainFile::parse_chain_line(
    std::string line, std::string &source, std::string &target,
    int &source_offset, int &target_offset, bool &same_strand
) {
    // Split `line` using space as deliminater.
    std::regex s_re("\\s+"); // space
    std::vector<std::string> vec(
        std::sregex_token_iterator(
            line.begin(), line.end(), s_re, -1),
        std::sregex_token_iterator());
    
    if (vec[0] == "chain"){
        if (vec.size() != 13){
            std::cerr << "Error during parsing chain" << line << "\n";
            std::cerr << "Invalid chain header\n";
            exit(1);
        }
        source = vec[2];
        target = vec[7];
        int source_length = stoi(vec[3]);
        same_strand = (vec[4] == vec[9]);
        try {
            source_offset = stoi(vec[5]);
            target_offset = stoi(vec[10]);
        }
        catch(...) {
            std::cerr << "Error during parsing chain" << line << "\n";
            exit(1);
        }
        this->init_bitvectors(source, source_length);
    } else if (vec.size() == 3) {
        int s_int_start = source_offset;
        int s_int_end = source_offset + stoi(vec[0]);
        int t_int_start = target_offset;
        int t_int_end = target_offset + stoi(vec[0]);
        // Set bitvectors
        if (s_int_start > 0)
            this->start_bv_map[source][s_int_start-1] = 1;
        this->end_bv_map[source][s_int_end] = 1;

        if (this->verbose > 1)
            fprintf(stderr, "source (%d-%d), target (%d-%d)\n", s_int_start, s_int_end, t_int_start, t_int_end);
        Interval intvl(target, s_int_start, s_int_end, t_int_start-s_int_start, same_strand);
        if (this->verbose > 1)
            intvl.debug_print_interval();

        std::unordered_map<std::string, std::vector<chain::Interval>>::const_iterator find_start = this->interval_map.find(source);
        if (find_start == this->interval_map.end()) {
            std::vector<chain::Interval> new_intervals;
            this->interval_map[source] = new_intervals;
        }
        this->interval_map[source].push_back(intvl);

        // Update offsets contributed by the variant
        source_offset = s_int_end + stoi(vec[1]);
        target_offset = t_int_end + stoi(vec[2]);
    } else {
        source = "";
        target = "";
        source_offset = 0;
        target_offset = 0;
        same_strand = true;
    }

    if (this->verbose > 1) {
        for (auto&& s: vec)
            std::cerr << s << " ";
        std::cerr << "\n";
    }
}


void ChainFile::show_interval_info(std::string contig, int pos) {
    int rank = this->get_start_rank(contig, pos);
    this->interval_map[contig][rank].debug_print_interval();
}

int ChainFile::get_lifted_pos(std::string contig, int pos) {
    int rank = this->get_start_rank(contig, pos);
    std::cerr << this->interval_map[contig][rank].target << "\n";
    int lifted_pos = pos + this->interval_map[contig][rank].offset;
    std::cerr << lifted_pos << "\n";
    return lifted_pos;
}

ChainFile* chain::chain_open(std::string fname, int verbose) {
    ChainFile *file = new ChainFile(fname, verbose);

    if (verbose > 0) {
        std::cerr << "Check if there are overlapping intervals...";
        std::string is_overlap = (file->interval_map_sanity_check())? " pass" : " failed";
        std::cerr << is_overlap << "\n";
    }

    // TODO
    std::cerr << "TEST_RANK\n";
    std::cerr << file->get_start_rank("chrM", 15000) << "\n";
    // std::cerr << file->get_start_rank("chrY", 6436563) << "\n";
    // std::cerr << file->get_start_rank("chrY", 9741965) << "\n";
    // std::cerr << file->get_start_rank("chrZ", 6436563) << "\n";

    std::cerr << "TEST1\n";
    int test_pos1 = 16000000;
    int test_rank1 = file->get_start_rank("chr1", test_pos1);
    std::cerr << "rank(" << test_pos1 << ")=" << test_rank1 << "\n";
    file->show_interval_info("chr1", test_pos1);
    file->get_lifted_pos("chr1", test_pos1);

    std::cerr << "TEST2\n";
    int test_pos2 = 674047;
    int test_rank2 = file->get_start_rank("chr1", test_pos2);
    std::cerr << "rank(" << test_pos2 << ")=" << test_rank2 << "\n";
    file->show_interval_info("chr1", test_pos2);
    file->get_lifted_pos("chr1", test_pos2);

    std::cerr << "TEST_RANK_END\n";

    return file;
}
