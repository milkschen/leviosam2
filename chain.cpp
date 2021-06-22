#include <iostream>
#include "chain.hpp"

using namespace chain;

size_t Interval::serialize(std::ofstream& out) const {
    size_t size = 0;
    size_t str_size = target.size();
    out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
    out.write(reinterpret_cast<const char*>(target.data()), str_size);
    size += sizeof(str_size) + str_size;

    out.write(reinterpret_cast<const char*>(&offset), sizeof(offset));
    out.write(reinterpret_cast<const char*>(&source_start), sizeof(source_start));
    out.write(reinterpret_cast<const char*>(&source_end), sizeof(source_end));
    out.write(reinterpret_cast<const char*>(&same_strand), sizeof(same_strand));
    size += sizeof(offset) + sizeof(source_start) + sizeof(source_end) + sizeof(same_strand);
    return size;
}

// Load from stream
void Interval::load(std::istream& in) {
    // target
    size_t str_size;
    in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
    std::vector<char> buf(str_size);
    in.read(reinterpret_cast<char*>(buf.data()), str_size);
    std::string target(buf.begin(), buf.end());
    int offset;
    in.read(reinterpret_cast<char*>(&offset), sizeof(offset));
    int source_start;
    in.read(reinterpret_cast<char*>(&source_start), sizeof(source_start));
    int source_end;
    in.read(reinterpret_cast<char*>(&source_end), sizeof(source_end));
    bool same_strand;
    in.read(reinterpret_cast<char*>(&same_strand), sizeof(same_strand));
    this->target = target;
    this->offset = offset;
    this->source_start = source_start;
    this->source_end = source_end;
    this->same_strand = same_strand;
    // std::cerr << "target=" << target << "\noffset=" << offset <<
    //     "\n(start, end, strand)=(" << source_start << "," << source_end << "," << same_strand << ")\n";
}

ChainMap::ChainMap(std::string fname, int verbose) {
    this->verbose = verbose;
    std::ifstream chain_f(fname);
    std::string line;
    std::unordered_map<std::string, sdsl::bit_vector> start_bv_map;
    std::unordered_map<std::string, sdsl::bit_vector> end_bv_map;
    if (chain_f.is_open()) {
        std::string source;
        std::string target;
        int source_offset = 0;
        int target_offset = 0;
        bool current_ss = true;
        while (getline(chain_f, line)) {
            this->parse_chain_line(
                line, source, target, source_offset, 
                target_offset, current_ss,
                start_bv_map, end_bv_map);
        }
        chain_f.close();
    }
    // biv_vector to sd_vector
    for (auto& it: start_bv_map) {
        auto itr = this->start_map.find(it.first);
        if (itr == this->start_map.end()) {
            this->start_map[it.first] = sdsl::sd_vector<>(it.second);
        }
    }
    for (auto& it: end_bv_map) {
        auto itr = this->end_map.find(it.first);
        if (itr == this->end_map.end()) {
            this->end_map[it.first] = sdsl::sd_vector<>(it.second);
        }
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
void ChainMap::init_bitvectors(
    std::string source, int source_length,
    std::unordered_map<std::string, sdsl::bit_vector> &start_bv_map,
    std::unordered_map<std::string, sdsl::bit_vector> &end_bv_map
) {
    auto itr = start_bv_map.find(source);
    if (itr == start_bv_map.end()) {
        sdsl::bit_vector new_start_bv(source_length), new_end_bv(source_length);
        start_bv_map[source] = new_start_bv;
        end_bv_map[source] = new_end_bv;
    }
}


void ChainMap::sort_interval_map() {
    for (auto &itr: this->interval_map) {
        this->sort_intervals(itr.first);
    }
}


/* Sort the intervals in a ChainMap in ascending order. */
void ChainMap::sort_intervals(std::string contig) {
    std::sort(
        this->interval_map[contig].begin(), this->interval_map[contig].end(),
        [](const Interval& lhs, const Interval& rhs) {
            return lhs.source_start < rhs.source_start;}
    );
}


void ChainMap::debug_print_interval_map() {
    for (auto &intvl_map: this->interval_map) {
        this->debug_print_intervals(intvl_map.first);
    }
}


/* Print the intervals in a ChainMap. */
void ChainMap::debug_print_intervals(std::string contig) {
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
bool ChainMap::interval_map_sanity_check() {
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
int ChainMap::get_start_rank(std::string contig, int pos) {
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
int ChainMap::get_end_rank(std::string contig, int pos) {
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
void ChainMap::init_rs() {
    for (auto& it: this->start_map) {
        std::string contig = it.first;
        std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_start = this->start_map.find(contig);
        if (find_start == this->start_map.end()) {
            // this->start_map[contig] = sdsl::sd_vector<>(it.second);
            this->start_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(this->start_rs1_map[contig], &this->start_map[contig]);
    }
    for (auto& it: this->end_map) {
        std::string contig = it.first;
        std::unordered_map<std::string, sdsl::sd_vector<>>::const_iterator find_start = this->end_map.find(contig);
        if (find_start == this->end_map.end()) {
            // this->end_map[contig] = sdsl::sd_vector<>(it.second);
            this->end_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(this->end_rs1_map[contig], &this->end_map[contig]);
    }
}

/* Parse a chain line and update the chain object.
 */
void ChainMap::parse_chain_line(
    std::string line, std::string &source, std::string &target,
    int &source_offset, int &target_offset, bool &same_strand,
    std::unordered_map<std::string, sdsl::bit_vector> &start_bv_map,
    std::unordered_map<std::string, sdsl::bit_vector> &end_bv_map
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
        this->init_bitvectors(source, source_length, start_bv_map, end_bv_map);
    } else if (vec.size() == 3) {
        int s_int_start = source_offset;
        int s_int_end = source_offset + stoi(vec[0]);
        int t_int_start = target_offset;
        int t_int_end = target_offset + stoi(vec[0]);
        // Set bitvectors
        if (s_int_start > 0)
            start_bv_map[source][s_int_start-1] = 1;
            // this->start_bv_map[source][s_int_start-1] = 1;
        end_bv_map[source][s_int_end] = 1;
        // this->end_bv_map[source][s_int_end] = 1;

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


void ChainMap::show_interval_info(std::string contig, int pos) {
    int rank = this->get_start_rank(contig, pos);
    this->interval_map[contig][rank].debug_print_interval();
}


std::string ChainMap::lift_contig(std::string contig, size_t pos) {
    int intvl_idx = this->get_start_rank(contig, pos);
    if (intvl_idx == -1)
        return "*";
    else
        return this->interval_map[contig][intvl_idx].target;
}

std::string ChainMap::lift_contig(
    std::string contig, int start_intvl_idx, int end_intvl_idx) {
    if (start_intvl_idx == -1)
        return "*";
    else
        return this->interval_map[contig][start_intvl_idx].target;
}

void ChainMap::lift_cigar(const std::string& contig, bam1_t* aln) {
    // TODO
}

size_t ChainMap::lift_pos(
    std::string contig, size_t pos) {
    int intvl_idx = this->get_start_rank(contig, pos);
    if (intvl_idx == -1)
        return pos;
    else
        return pos + this->interval_map[contig][intvl_idx].offset;
}

size_t ChainMap::lift_pos(
    std::string contig, size_t pos,
    int start_intvl_idx, int end_intvl_idx) {
    if (start_intvl_idx == -1)
        return pos;
    else
        return pos + this->interval_map[contig][start_intvl_idx].offset;
}

// Return liftability (bool) and update start/end-interval indexes.
bool ChainMap::is_liftable(
    std::string contig, size_t pos,
    int &start_intvl_idx, int &end_intvl_idx
) {
    start_intvl_idx = this->get_start_rank(contig, pos);
    end_intvl_idx = this->get_end_rank(contig, pos);
    if (start_intvl_idx == -1 || end_intvl_idx == -1)
        return false;
    if (end_intvl_idx - start_intvl_idx == 1)
        return true;
    return false;
}

bool ChainMap::lift_segment(
    bam1_t* aln, bam_hdr_t* hdr,
    bool is_first_seg, std::string &dest_contig,
    int &start_intvl_idx, int &end_intvl_idx
) {
    bam1_core_t* c = &(aln->core);

    if (is_first_seg && (c->flag & BAM_FUNMAP))
        return false;
    else if (!is_first_seg && (c->flag & BAM_FMUNMAP))
        return false;

    std::string source_contig;
    if (is_first_seg)
        source_contig = hdr->target_name[c->tid];
    else
        source_contig = hdr->target_name[c->mtid];

    bool is_liftable;
    // Check liftability.
    // If not liftable (un-liftable or unmapped), 
    // update start/end_intvl_idx and return false 
    if (is_first_seg) {
        start_intvl_idx = this->get_start_rank(source_contig, c->pos);
        end_intvl_idx = this->get_end_rank(source_contig, c->pos);
    } else {
        start_intvl_idx = this->get_start_rank(source_contig, c->mpos);
        end_intvl_idx = this->get_end_rank(source_contig, c->mpos);
    }
    if ((start_intvl_idx == -1) || (end_intvl_idx == -1)) {
        return false;
    ////// #TODO } else if (start_intvl_idx - end_intvl_idx != 1) {
    ////// #TODO     return false;
    } else {
        is_liftable = true;
    }

    // Lift contig
    if (start_intvl_idx == -1)
        dest_contig = "*";
    else
        dest_contig = this->interval_map[source_contig][start_intvl_idx].target;
    if (is_first_seg)
        c->tid = sam_hdr_name2tid(hdr, dest_contig.c_str());
    else
        c->mtid = sam_hdr_name2tid(hdr, dest_contig.c_str());

    // Lift pos
    size_t pos;
    if (is_first_seg) {
        if (start_intvl_idx != -1)
            c->pos = c->pos +
                this->interval_map[source_contig][start_intvl_idx].offset;
    } else {
        if (start_intvl_idx != -1)
            c->mpos = c->mpos + 
                this->interval_map[source_contig][start_intvl_idx].offset;
    }

    // Lift cigar #TODO
    // this->lift_cigar(
    //     source_contig, aln, start_intvl_idx, end_intvl_idx);

    return is_liftable;
}

void ChainMap::lift_aln(
    bam1_t* aln,
    bam_hdr_t* hdr,
    std::string &dest_contig
) {
    bam1_core_t* c = &(aln->core);
    std::string source_contig;
    size_t lift_status;

    int start_intvl_idx, end_intvl_idx;
    size_t pos = c->pos;
    size_t mpos = c->mpos;
    bool is_r1_liftable = this->lift_segment(
        aln, hdr, true, dest_contig, start_intvl_idx, end_intvl_idx);

    std::string lo;
    // Paired
    if (c->flag & BAM_FPAIRED) {
        int mate_start_intvl_idx, mate_end_intvl_idx;
        bool is_r2_liftable = this->lift_segment(
            aln, hdr, false, dest_contig, start_intvl_idx, end_intvl_idx);

        // R1 unmapped
        if (c->flag & BAM_FUNMAP) {
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UM_R2_UM;
                // Neither lifted - do nothing
                lo = "UM_UM";
            } else if (is_r2_liftable) {
                lift_status = LIFT_R1_UM_R2_L;
                // TODO: copy R2 to R1
                c->pos = c->mpos;
                c->tid = c->mtid;
                lo = "UM_L";
            } else {
                lift_status = LIFT_R1_UM_R2_UL;
                // Neither lifted - do nothing
                lo = "UM_UL";
            }
        // R1 mapped, liftable
        } else if (is_r1_liftable) {
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_L_R2_UM;
                // TODO: Copy r1 to r2
                c->mpos = c->pos;
                c->mtid = c->tid;
                lo = "L_UM";
            } else if (is_r2_liftable) {
                lift_status = LIFT_R1_L_R2_L;
                lo = "L_L";
            } else {
                lift_status = LIFT_R1_L_R2_UL;
                // TODO: Copy r1 to r2
                c->mpos = c->pos;
                c->mtid = c->tid;
                lo = "L_UL";
            }
        // R1 mapped, un-liftable
        } else {
            // c->qual = 255;
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UL_R2_UM;
                // Neither lifted - do nothing
                lo = "UL_UM";
            } else if (is_r2_liftable) {
                lift_status = LIFT_R1_UL_R2_L;
                // #TODO Copy R2 to R1 ??
                c->pos = c->mpos;
                c->tid = c->mtid;
                lo = "UL_L";
            } else {
                lift_status = LIFT_R1_UL_R2_UL;
                // Neither lifted - do nothing
                lo = "UL_UL";
            }
        }
        // TODO consider STRAND
        c->isize = 
            (c->isize == 0)? 0 : 
            (c->tid != c->mtid)? 0 : 
                c->isize + (mpos - c->mpos) - (pos - c->pos);
    // Unpaired
    } else {
        if (c->flag & BAM_FUNMAP) {
            lift_status = LIFT_R_UM;
            lo = "UM";
        } else if (is_r1_liftable) {
            lift_status = LIFT_R_L;
            lo = "L";
        } else {
            lift_status = LIFT_R_UL;
            lo = "UL";
        }
    }
    bam_aux_append(
        aln, "LO", 'Z', lo.length() + 1,
        reinterpret_cast<uint8_t*>(const_cast<char*>(lo.c_str())));
}

// Save to stream
size_t ChainMap::serialize(std::ofstream& out) {
    size_t size = 0;
    // interval_map
    size_t nelems = interval_map.size();
    std::cerr << "serializing interval maps...\n";
    size += sizeof(nelems);
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    for (auto& x: this->interval_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        size += sizeof(str_size) + str_size;
        size_t vec_size = x.second.size();
        out.write(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
        size += sizeof(vec_size);
        for (auto &intvl: x.second) {
            size += intvl.serialize(out);
        }
    }
    // start_map
    std::cerr << "serializing start_maps...\n";
    nelems = this->start_map.size();
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    size += sizeof(nelems);
    for (auto& x: this->start_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        size += sizeof(str_size) + str_size;
        size += x.second.serialize(out);
    }
    // end_map
    std::cerr << "serializing end_maps...\n";
    nelems = this->end_map.size();
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    size += sizeof(nelems);
    for (auto& x: this->end_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        size += sizeof(str_size) + str_size;
        size += x.second.serialize(out);
    }
    return size;
}

ChainMap::ChainMap(std::ifstream& in) {
    this->load(in);
    // debug messages
    // this->debug_print_interval_map();
}

// loads from stream
void ChainMap::load(std::ifstream& in) {
    size_t map_size;
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        size_t vec_size;
        in.read(reinterpret_cast<char*>(&vec_size), sizeof(vec_size));
        for (auto j = 0; j < vec_size; ++j) {
            this->interval_map[key].push_back(Interval(in));
        }
    }
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        sdsl::sd_vector<> sdv;
        sdv.load(in);
        this->start_map[key] = sdv;
    }
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        sdsl::sd_vector<> sdv;
        sdv.load(in);
        this->end_map[key] = sdv;
    }
    this->init_rs();
}

