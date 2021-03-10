#include <iostream>
#include "chain.hpp"

using namespace chain;

bool DEBUG = false;

ChainFile::ChainFile(std::string fname) {
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
    if (DEBUG) {
        this->debug_print_interval_map();
        std::cerr << "\n\nsort intervals\n\n";
    }
    this->sort_interval_map();
    if (DEBUG) {
        this->debug_print_interval_map();
    }
}

// Create start and end bitvectors when see a new `source`.
// Do nothing if `source` has been seen.
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

// Sort the intervals in a ChainFile in ascending order.
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

// Print the intervals in a ChainFile.
void ChainFile::debug_print_intervals(std::string contig) {
    for (auto &intvl: this->interval_map[contig]) {
        intvl.debug_print_interval();
    }
    std::cerr << "\n";
}

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
        this->start_bv_map[source][s_int_start] = 1;
        this->end_bv_map[source][s_int_end] = 1;

        if (DEBUG)
            fprintf(stderr, "source (%d-%d), target (%d-%d)\n", s_int_start, s_int_end, t_int_start, t_int_end);
        Interval intvl(target, s_int_start, s_int_end, t_int_start-s_int_start, same_strand);
        // Interval intvl(source, target, s_int_start, s_int_end, t_int_start-s_int_start, same_strand);
        if (DEBUG)
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

    // Check if there are overlapping intervals
    // std::unordered_map<std::string, std::vector<chain::Interval>>::const_iterator imap_start = this->interval_map.begin();
    // std::unordered_map<std::string, std::vector<chain::Interval>>::iterator itr;
    // for (itr = this->interval_map.begin(); itr != this->interval_map.end(); itr++) {

    // for (auto& itr : this->interval_map){
    //         std::vector<chain::Interval> v = this->interval_map[itr.first];
    //         for (int i = 0; i < v.size()-2; i++) {
    //             v[i].debug_print_interval();
    //             if (v[i].source_start == v[i+1].source_start) {
    //                 std::cerr << "error: " << itr.first << "\n";
    //             }
    //         }
    // }
    
    std::vector<chain::Interval> v = this->interval_map["chr2"];
    for (int i = 0; i < 10; i++) {
        v[i].debug_print_interval();
    }

    if (DEBUG) {
        for (auto&& s: vec)
            std::cerr << s << " ";
        std::cerr << "\n";
    }
}

ChainFile* chain::chain_open(std::string fname) {
    ChainFile *file = new ChainFile(fname);
    
    // TODO
    std::cerr << "TEST_RANK\n";
    std::cerr << file->get_start_rank("chrM", 16000) << "\n";
    std::cerr << file->get_start_rank("chr1", 16000000) << "\n";
    // std::cerr << file->get_start_rank("chrY", 6436563) << "\n";
    std::cerr << file->get_start_rank("chrY", 9741965) << "\n";
    std::cerr << "TEST_RANK_END\n";

    return file;
}
