#include <cmath>
#include <iostream>
#include <queue>
#include "chain.hpp"


namespace chain {
Interval::Interval() {
    target = "";
    offset = 0;
    source_start = 0;
    source_end = 0;
    strand = true;
}

Interval::Interval(
    std::string t, int so, int se, int o, bool ss
) {
    target = t;
    offset = o;
    source_start = so;
    source_end = se;
    strand = ss;
}

Interval::Interval(std::ifstream& in) {
    this->load(in);
}

void Interval::debug_print_interval() {
    std::string ss = (strand)? "+" : "-";
    std::cerr << "[" << source_start << ":" << source_end << ")->" << target << " (" << ss << "); offset = " << offset << "\n";
}

size_t Interval::serialize(std::ofstream& out) const {
    size_t size = 0;
    size_t str_size = target.size();
    out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
    out.write(reinterpret_cast<const char*>(target.data()), str_size);
    size += sizeof(str_size) + str_size;

    out.write(reinterpret_cast<const char*>(&offset), sizeof(offset));
    out.write(reinterpret_cast<const char*>(&source_start), sizeof(source_start));
    out.write(reinterpret_cast<const char*>(&source_end), sizeof(source_end));
    out.write(reinterpret_cast<const char*>(&strand), sizeof(strand));
    size += sizeof(offset) + sizeof(source_start) + sizeof(source_end) + sizeof(strand);
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
    bool strand;
    in.read(reinterpret_cast<char*>(&strand), sizeof(strand));
    this->target = target;
    this->offset = offset;
    this->source_start = source_start;
    this->source_end = source_end;
    this->strand = strand;
    // std::cerr << "target=" << target << "\noffset=" << offset <<
    //     "\n(start, end, strand)=(" << source_start << "," << source_end << "," << strand << ")\n";
}


ChainMap::ChainMap(std::string fname, int verbose) {
    this->verbose = verbose;
    std::ifstream chain_f(fname);
    std::string line;
    BitVectorMap start_bv_map;
    BitVectorMap end_bv_map;
    if (chain_f.is_open()) {
        std::string source;
        std::string target;
        int source_offset = 0;
        int target_offset = 0;
        int source_len = 0;
        bool current_ss = true;
        while (getline(chain_f, line)) {
            this->parse_chain_line(
                line, source, target, source_len, source_offset, 
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
    if (this->verbose > 2) {
        this->debug_print_interval_map();
        std::cerr << "\n\nsort intervals\n\n";
    }
    this->sort_interval_map();
    if (this->verbose > 2) {
        this->debug_print_interval_map();
    }
    // TEMP
    this->interval_map_sanity_check();
}


/* Create start and end bitvectors when see a new `source`.
 * Do nothing if `source` has been seen.
 */
void ChainMap::init_bitvectors(
    std::string source, int source_len,
    BitVectorMap &start_bv_map, BitVectorMap &end_bv_map
) {
    auto itr = start_bv_map.find(source);
    if (itr == start_bv_map.end()) {
        sdsl::bit_vector new_start_bv(source_len), new_end_bv(source_len);
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
    std::cerr << "Interval_map sanity check: passed (no overlaps)\n";
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
    SdVectorMap::const_iterator find_start = this->start_map.find(contig);
    if (find_start == this->start_map.end())
        return -1;
    if (pos >= start_rs1_map[contig].size())
        pos = start_rs1_map[contig].size() - 1;
    return start_rs1_map[contig](pos);
}

/* Get the rank in the end bitvector at contig[pos]
 * If pos > len(contig): return rank(contig[-1])
 * Returns:
 *   non-negative int: rank
 *   -1: if pos < 0 or contig does not exist
 */
int ChainMap::get_end_rank(std::string contig, int pos) {
    if (pos < 0)
        return -1;
    SdVectorMap::const_iterator find_end = this->end_map.find(contig);
    if (find_end == this->end_map.end())
        return -1;
    if (pos >= end_rs1_map[contig].size())
        pos = end_rs1_map[contig].size() - 1;
    return end_rs1_map[contig](pos);
}

/* Init rank support for all start/end bitvectors. */
void ChainMap::init_rs() {
    for (auto& it: this->start_map) {
        std::string contig = it.first;
        auto itr = this->start_map.find(contig);
        if (itr == this->start_map.end()) {
            this->start_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(this->start_rs1_map[contig], &this->start_map[contig]);
    }
    for (auto& it: this->end_map) {
        std::string contig = it.first;
        auto itr = this->end_map.find(contig);
        if (itr == this->end_map.end()) {
            this->end_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(this->end_rs1_map[contig], &this->end_map[contig]);
    }
}

/* Parse a chain line and update the chain object.
 */
void ChainMap::parse_chain_line(
    std::string line, std::string &source, std::string &target,
    int &source_len, int &source_offset, int &target_offset, bool &strand,
    BitVectorMap &start_bv_map, BitVectorMap &end_bv_map
) {
    // Split `line` using space as deliminater.
    std::regex s_re("\\s+"); // space
    std::vector<std::string> vec(
        std::sregex_token_iterator(
            line.begin(), line.end(), s_re, -1),
        std::sregex_token_iterator());
    
    int target_len;
    if (vec[0] == "chain"){
        if (vec.size() != 13){
            std::cerr << "Error during parsing chain" << line << "\n";
            std::cerr << "Invalid chain header\n";
            exit(1);
        }
        source = vec[2];
        target = vec[7];
        source_len = stoi(vec[3]);
        target_len = stoi(vec[8]);
        strand = (vec[4] == vec[9]);
        // std::cerr << source << "->" << target << " " << target_len << " " << strand << "\n";
        try {
            source_offset = stoi(vec[5]);
            if (strand) {
                target_offset = stoi(vec[10]);
            } else {
                // target_offset = target_len - stoi(vec[11]);
                target_offset = target_len - stoi(vec[10]);
            }
        }
        catch(...) {
            std::cerr << "Error during parsing chain:\n" << line << "\n";
            exit(1);
        }
        this->init_bitvectors(source, source_len, start_bv_map, end_bv_map);
    } else if (vec[0] != "") {
        int s_int_start = source_offset;
        int t_int_start = target_offset;
        int s_int_end = source_offset + stoi(vec[0]);
        int t_int_end;
        if (strand) {
            t_int_end = target_offset + stoi(vec[0]);
        } else {
            t_int_end = target_offset - stoi(vec[0]);
        }
        // Set bitvectors
        if (s_int_start > 0)
            start_bv_map[source][s_int_start-1] = 1;
        else if (s_int_start == 0)
            start_bv_map[source][0] = 1;
        if (s_int_end >= source_len)
            end_bv_map[source][source_len-1] = 1;
        else
            // end_bv_map[source][s_int_end-1] = 1;
            end_bv_map[source][s_int_end] = 1;

        if (this->verbose > 1)
            fprintf(stderr, "source (%d-%d), target (%d-%d)\n",
                    s_int_start, s_int_end, t_int_start, t_int_end);
        int offset;
        if (strand) {
            offset = t_int_start - s_int_start;
            // Interval intvl(target, s_int_start, s_int_end,
            //                t_int_start-s_int_start, strand);
        } else {
            offset = t_int_end - s_int_start;
        }
        Interval intvl(target, s_int_start, s_int_end,
                       offset, strand);
        if (this->verbose > 1)
            intvl.debug_print_interval();

        IntervalMap::const_iterator find_start =
            this->interval_map.find(source);
        if (find_start == this->interval_map.end()) {
            std::vector<chain::Interval> new_intervals;
            this->interval_map[source] = new_intervals;
        }
        this->interval_map[source].push_back(intvl);

        if (vec.size() == 3) {
            // Update offsets contributed by the variant
            source_offset = s_int_end + stoi(vec[1]);
            if (strand) {
                target_offset = t_int_end + stoi(vec[2]);
            } else {
                target_offset = t_int_end - stoi(vec[2]);
            }
        }
    } else {
        source = "";
        target = "";
        source_offset = 0;
        target_offset = 0;
        strand = true;
    }

    if (this->verbose > 1) {
        for (auto&& s: vec)
            std::cerr << s << " ";
        std::cerr << "\n";
    }
}


/* Update SAM flag to set a record as an unmapped alignment
 *   - Clear forward/reverse status.
 *   - If paired, changed to improper paired. 
 */
void ChainMap::update_flag_unmap(bam1_core_t* c, const bool is_first_seg) {
    if (is_first_seg) {
        c->flag |= BAM_FUNMAP;
        c->flag &= ~BAM_FPROPER_PAIR;
        c->flag &= ~BAM_FREVERSE;
    } else {
        c->flag |= BAM_FMUNMAP;
        c->flag &= ~BAM_FPROPER_PAIR;
        c->flag &= ~BAM_FMREVERSE;
    }
}

/* Lift CIGAR
 * This version of `lift_cigar()` can be called when no interval indexes are
 * available, but is slower because of additional calls to query for interval indexes.
 */
void ChainMap::lift_cigar(const std::string& contig, bam1_t* aln) {
    bam1_core_t* c = &(aln->core);
    auto pos_end = c->pos + bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));
    // Note that intvl_idx = rank - 1
    auto start_intvl_idx = this->get_start_rank(contig, c->pos) - 1;
    auto end_intvl_idx = this->get_end_rank(contig, c->pos) - 1;
    auto pend_start_intvl_idx = this->get_start_rank(contig, pos_end) - 1;

    auto next_intvl = this->interval_map[contig][start_intvl_idx+1];
    int num_sclip_start = 0;
    if ((start_intvl_idx == end_intvl_idx) &&
        (start_intvl_idx < this->interval_map[contig].size() - 1))
        num_sclip_start = std::abs(c->pos - next_intvl.source_start);
    this->lift_cigar(
        contig, aln,
        start_intvl_idx, pend_start_intvl_idx, num_sclip_start);
}

/* Lift CIGAR
 * A fast version of `lift_cigar()`, which is used when auxiliary information 
 * including query start/end interval ranks and #clipped bases is available.
 */
void ChainMap::lift_cigar(
    const std::string& contig, bam1_t* aln,
    int start_intvl_idx, int pend_start_intvl_idx, int num_clipped) {

    // If POS is inside `current_intvl`, `num_clipped` <= 0
    if ((num_clipped <= 0) && (start_intvl_idx == pend_start_intvl_idx)) {
        return;
    }
    // auto current_intvl = this->interval_map[contig][start_intvl_idx];

    bam1_core_t* c = &(aln->core);
    uint32_t* cigar = bam_get_cigar(aln);
    std::vector<uint32_t> new_cigar;
    lift_cigar_core(new_cigar, contig, aln, start_intvl_idx, pend_start_intvl_idx, num_clipped);

    // Adjust data position when n_cigar is changed by levioSAM. n_cigar could either be
    // increased or decreased.
    //
    // Adapted from samtools/sam.c
    // https://github.com/samtools/htslib/blob/2264113e5df1946210828e45d29c605915bd3733/sam.c#L515
    if (aln->core.n_cigar != new_cigar.size()){
    // if (aln->core.n_cigar != cigar_op_len.size()){
        auto cigar_st = (uint8_t*)bam_get_cigar(aln) - aln->data;
        auto fake_bytes = aln->core.n_cigar * 4;
        aln->core.n_cigar = (uint32_t)new_cigar.size();
        // aln->core.n_cigar = (uint32_t)cigar_op_len.size();
        auto n_cigar4 = aln->core.n_cigar * 4;
        auto orig_len = aln->l_data;
        if (n_cigar4 > fake_bytes){
            // Check if we need to update `aln->m_data`.
            //
            // Adapted from htslib/sam_internal.h
            // https://github.com/samtools/htslib/blob/31f0a76d338c9bf3a6893b71dd437ef5fcaaea0e/sam_internal.h#L48
            auto new_m_data = (size_t) aln->l_data + n_cigar4 - fake_bytes;
            kroundup32(new_m_data);
            if (new_m_data > aln->m_data){
                auto new_data = static_cast<uint8_t *>(realloc(aln->data, new_m_data));
                if (!new_data){
                    std::cerr << "[Error] Failed to expand a bam1_t struct for " << bam_get_qname(aln) << "\n";
                    std::cerr << "This is likely due to out of memory\n";
                    exit(1);
                }
                aln->data = new_data;
                aln->m_data = new_m_data;
                cigar = bam_get_cigar(aln);
            }
        }
        // Update memory usage of data.
        aln->l_data = aln->l_data - fake_bytes + n_cigar4;
        // Move data to the correct place.
        memmove(aln->data + cigar_st + n_cigar4,
                aln->data + cigar_st + fake_bytes,
                orig_len - (cigar_st + fake_bytes));
        // If new n_cigar is greater, copy the real CIGAR to the right place.
        // Skipped this if new n_cigar is smaller than the original value.
        if (n_cigar4 > fake_bytes){
            memcpy(aln->data + cigar_st,
                   aln->data + (n_cigar4 - fake_bytes) + 8,
                   n_cigar4);
        }
        aln->core.n_cigar = new_cigar.size();
        // aln->core.n_cigar = cigar_op_len.size();
    }
    for (int i = 0; i < aln->core.n_cigar; i++){
        *(cigar + i) = new_cigar[i];
    }
    if (this->verbose > 3) {
        std::cerr << "* new: ";
        debug_print_cigar(aln);
    }
}

void ChainMap::lift_cigar_core(
    std::vector<uint32_t>& new_cigar, const std::string& contig,
    bam1_t* aln, int start_intvl_idx, int pend_start_intvl_idx, int num_clipped) {
    uint32_t* cigar = bam_get_cigar(aln);
    bam1_core_t* c = &(aln->core);

    bool TMP_DEBUG = (this->verbose > 3);
    // DEBUG
    if (TMP_DEBUG) {
        std::cerr << "\n" << bam_get_qname(aln) << "\n";
        std::cerr << "  Num_clipped = " << num_clipped << "\n";
        std::cerr << "    start " << start_intvl_idx << ", pend " << pend_start_intvl_idx << "\n";
        std::cerr << "* old: ";
        for (int i = 0; i < aln->core.n_cigar; i++){
            auto cigar_op_len = bam_cigar_oplen(cigar[i]);
            auto cigar_op = bam_cigar_op(cigar[i]);
            std::cerr << cigar_op_len << bam_cigar_opchr(cigar_op);
        }
        std::cerr << "\n";
        std::cerr << "pos: " << c->pos << "\n";
    }
    // END_DEBUG

    std::queue<std::tuple<int, int>> break_points;
    for (auto i = start_intvl_idx; i <= pend_start_intvl_idx; i++) {
        int bp = this->interval_map[contig][i].source_end - c->pos;
        int diff = this->interval_map[contig][i+1].offset -
                   this->interval_map[contig][i].offset;
        if (bp > 0) {
            break_points.push(std::make_tuple(bp, diff));
            // DEBUG
            if (TMP_DEBUG) {
                this->interval_map[contig][i].debug_print_interval();
                std::cerr << "bp=" << bp << ", diff=" << diff << "\n";
            }
            // END_DEBUG
        }
        // break_points.push(std::make_tuple(bp, diff));
    }

    // Number of bases that need to be clipped in next iterations
    int tmp_clipped = 0;
    num_clipped = (num_clipped < 0)? 0 : num_clipped;
    int query_offset = 0;
    int tmp_gap = 0;
    for (auto i = 0; i < aln->core.n_cigar; i++){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        if (TMP_DEBUG) {
            std::cerr << "  OP=" << cigar_op << ", OP_LEN=" << cigar_op_len << "\n";
        }
        if (num_clipped == 0) {
            // If within one interval, update CIGAR and jump to the next CIGAR operator
            if (start_intvl_idx == pend_start_intvl_idx) {
                new_cigar.push_back(bam_cigar_gen(cigar_op_len, cigar_op));
                continue;
            }

            if (bam_cigar_type(cigar_op) & 1) {
                // Resolve tmp gap bases
                if (tmp_gap > 0) {
                    if (tmp_gap < cigar_op_len) {
                        new_cigar.push_back(bam_cigar_gen(tmp_gap, BAM_CDEL));
                        cigar_op_len -= tmp_gap;
                        tmp_gap = 0;
                    } else {
                        new_cigar.push_back(bam_cigar_gen(cigar_op_len, BAM_CDEL));
                        tmp_gap -= cigar_op_len;
                        continue;
                    }
                }
                auto next_bp = std::get<0>(break_points.front());
                auto next_q_offset = query_offset + cigar_op_len;
                if (TMP_DEBUG) {
                    std::cerr << "next_bp =" << next_bp << ", next_q_offset =" << next_q_offset << "\n";
                    std::cerr << "original CIGAR = " << cigar_op_len << bam_cigar_opchr(cigar_op) << "\n";
                }
                // Have not reached the next breakpoint
                if (next_q_offset <= next_bp){
                    if (TMP_DEBUG) {
                        std::cerr << "have not reach next_bp: " << cigar_op_len << bam_cigar_opchr(cigar_op) << "\n";
                    }
                    new_cigar.push_back(bam_cigar_gen(cigar_op_len, cigar_op));
                // Split one CIGAR chunk into two parts and insert lift-over bases there
                } else {
                    auto curr_bp = next_bp;
                    auto num_ins = 0;
                    int second_half_len = cigar_op_len;
                    while (next_bp >= query_offset && next_q_offset > next_bp) {
                        curr_bp = next_bp;
                        auto first_half_len = next_bp - query_offset;
                        if (TMP_DEBUG) {
                            std::cerr << "first half: " << first_half_len << bam_cigar_opchr(cigar_op) << "\n";
                        }
                        if (first_half_len > 0)
                            new_cigar.push_back(bam_cigar_gen(first_half_len, cigar_op));
                        second_half_len -= first_half_len;
                        query_offset += first_half_len;
                        auto diff = std::get<1>(break_points.front());
                        if (diff > 0) { // D
                            new_cigar.push_back(bam_cigar_gen(diff, BAM_CDEL));
                            if (TMP_DEBUG) {
                                std::cerr << "update liftover CDEL: " << diff << bam_cigar_opchr(BAM_CDEL) << "\n";
                            }
                        } else { // I
                            // Adding `diff` makes the CIGAR longer than QUERY
                            if (query_offset - diff > c->l_qseq)
                                diff = query_offset - c->l_qseq;
                            new_cigar.push_back(bam_cigar_gen(-diff, BAM_CINS));
                            num_ins -= diff;
                            second_half_len += diff;
                            query_offset -= diff;
                            if (TMP_DEBUG) {
                                std::cerr << "update liftover CINS: " << -diff << bam_cigar_opchr(BAM_CINS) << "\n";
                            }
                        }

                        break_points.pop();
                        next_bp = std::get<0>(break_points.front());
                    }
                    if (TMP_DEBUG) {
                        std::cerr << "second_half_len=" << second_half_len << "\n";
                    }
                    if (second_half_len < 0) {
                        tmp_gap = -second_half_len;
                    } else if (second_half_len > 0) {
                        new_cigar.push_back(bam_cigar_gen(second_half_len, cigar_op));
                        if (TMP_DEBUG) {
                            std::cerr << second_half_len << bam_cigar_opchr(cigar_op) << "\n";
                        }
                    }
                }
                query_offset = next_q_offset;
                if (TMP_DEBUG) {
                    std::cerr << "query_offset -> " << query_offset << "\n";
                }
            } else {
                // TODO
                new_cigar.push_back(bam_cigar_gen(cigar_op_len, cigar_op));
            }
            // new_cigar.push_back(bam_cigar_gen(cigar_op_len, cigar_op));
            continue;
        }

        // Handle soft clipped bases
        // Only replace bases that consume the QUERY with the SOFT_CLIP ("S") operator
        if (bam_cigar_type(cigar_op) & 1) {
            // If there are bases need clipping, check if this OP has enough bases to clip:
            //   If so, add soft clippings ("S") and update the length of the current OP
            //   If not, cancel out OPLEN bases of clippings
            if (cigar_op_len >= num_clipped) {
                cigar_op_len -= num_clipped;
                new_cigar.push_back(bam_cigar_gen(num_clipped + tmp_clipped, BAM_CSOFT_CLIP));
                query_offset += num_clipped + tmp_clipped;
                num_clipped = 0;
                tmp_clipped = 0;
                if (cigar_op_len > 0) {
                    new_cigar.push_back(bam_cigar_gen(cigar_op_len, cigar_op));
                    query_offset += cigar_op_len;
                }
            } else {
                tmp_clipped += cigar_op_len;
                num_clipped -= cigar_op_len;
            }
        }
        // } else {
        //     if (cigar_op_len <= num_clipped) {
        //         tmp_clipped -= cigar_op_len;
        //         num_clipped += cigar_op_len;
        //     } else {
        //         new_cigar.push_back(
        //             bam_cigar_gen(cigar_op_len - num_clipped, cigar_op));
        //     }
        // }
    }
    if (tmp_gap > 0) {
        // TODO: trim bases from the end
    }
    if (TMP_DEBUG) {
        std::cerr << "  len(new_cigar) = " << new_cigar.size() << "\n";
    }
}

// Return false if unliftable
bool ChainMap::lift_segment(
    bam1_t* aln, bam_hdr_t* hdr,
    bool is_first_seg, std::string &dest_contig,
    int &start_intvl_idx, int &end_intvl_idx
) {
    const int allowed_num_indel = 10;
    bam1_core_t* c = &(aln->core);

    // If unmapped, the alignment is not liftable
    // When the mate is aligned, the info of the unmapped segment will be updated in lift_aln()
    if (is_first_seg && (c->flag & BAM_FUNMAP))
        return false;
    else if (!is_first_seg && (c->flag & BAM_FMUNMAP))
        return false;

    std::string source_contig = (is_first_seg)?
        hdr->target_name[c->tid] : hdr->target_name[c->mtid];

    if (is_first_seg) {
        // Note that intvl_idx = rank - 1
        start_intvl_idx = this->get_start_rank(source_contig, c->pos) - 1;
        end_intvl_idx = this->get_end_rank(source_contig, c->pos) - 1;
    } else {
        start_intvl_idx = this->get_start_rank(source_contig, c->mpos) - 1;
        end_intvl_idx = this->get_end_rank(source_contig, c->mpos) - 1;
    }

    /* Check liftability
     * If an aln cannot be mapped to an interval
     *   (a) the position is not valid in the intervalmap
     *   (b) the position is within the gap between two valid intervals,
     *   set it as unliftable.
     *
     * An exception is if the aln is close enough to the next interval (under scenario (b)),
     * we might choose to rescue it, with the bases in the gap clipped.
     */
    bool is_liftable;
    int num_sclip_start = 0;
    if ((start_intvl_idx <= -1) || (end_intvl_idx <= -1)) {
        update_flag_unmap(c, is_first_seg);
        return false;
    } else if (start_intvl_idx == end_intvl_idx) {
        /* If the starting position is not within any intervals in the ChainMap,
         * this segment is set to unmapped (BAM_FUNMAP or BAM_FMUNMAP).
         *
         * If the alignment is close enough with the next interval (`<allowed_num_indel`),
         * we still consier this alignment as liftable
         */
        if (start_intvl_idx < this->interval_map[source_contig].size() - 1) {
            auto next_intvl = this->interval_map[source_contig][start_intvl_idx+1];
            num_sclip_start = std::abs(c->pos - next_intvl.source_start);
            if (num_sclip_start >= allowed_num_indel) {
                update_flag_unmap(c, is_first_seg);
                return false;
            }
            // Advance start_intvl_idx
            start_intvl_idx += 1;
            is_liftable = true;
        } else {
            update_flag_unmap(c, is_first_seg);
            return false;
        }
    } else {
        is_liftable = true;
    }
    
    // Debug messages
    if (this->verbose > 1) {
        if (is_first_seg) {
            std::cerr << "First segment\n";
            std::cerr << "  Source: " << source_contig << ":" << c->pos << "\n";
        } else {
            std::cerr << "Second segment\n";
            std::cerr << "  Source: " << source_contig << ":" << c->mpos << "\n";
        }
        std::cerr << "  Queried interval indexes:\n";
        std::cerr << "    start: " << start_intvl_idx << "\n    ";
        this->interval_map[source_contig][start_intvl_idx].debug_print_interval();
        std::cerr << "    end  : " << end_intvl_idx << "\n";
    }


    /* Lift contig */
    auto current_intvl = this->interval_map[source_contig][start_intvl_idx];
    dest_contig = current_intvl.target;
    if (is_first_seg)
        c->tid = sam_hdr_name2tid(hdr, dest_contig.c_str());
    else
        c->mtid = sam_hdr_name2tid(hdr, dest_contig.c_str());

    // Lift pos
    auto offset = current_intvl.offset;
    auto strand = current_intvl.strand;
    auto pos_end = c->pos + bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));
    // Estimate ending pos of the mate
    auto mpos_end = c->mpos + c->l_qseq;
    // Note that intvl_idx = rank - 1
    auto pend_start_intvl_idx = (is_first_seg)? 
        this->get_start_rank(source_contig, pos_end) - 1 :
        this->get_start_rank(source_contig, mpos_end) - 1;


    /* If an alignment is overlapped with multiple intervals, check the gap size between the 
     * intervals. If the gap size is greater than `allowed_num_indel`, mark the alignment as
     * unliftable (return false).
     */
    if (start_intvl_idx != pend_start_intvl_idx) {
        for (auto j=start_intvl_idx; j < pend_start_intvl_idx - 1; j ++) {
            if (std::abs(this->interval_map[source_contig][j+1].offset -
                         this->interval_map[source_contig][j].offset) > allowed_num_indel) {
                update_flag_unmap(c, is_first_seg);
                return false;
            }
        }
        // auto next_intvl = this->interval_map[source_contig][pend_start_intvl_idx];
        // if (std::abs(next_intvl.offset - current_intvl.offset) > allowed_num_indel) {
        //     // is_liftable = false;
        //     update_flag_unmap(c, is_first_seg);
        //     return false;
        // }
    }

    // Debug messages
    if (this->verbose > 1) {
        std::cerr << "* Estimated ending position:\n";
        if (is_first_seg) {
            std::cerr << "  Source: " << source_contig << ":" << pos_end << "\n";
        } else {
            std::cerr << "  Source: " << source_contig << ":" << mpos_end << "\n";
        }
        std::cerr << "  Queried interval indexes:\n";
        std::cerr << "    start: " << pend_start_intvl_idx << "\n    ";
        this->interval_map[source_contig][pend_start_intvl_idx].debug_print_interval();
        std::cerr << "\n";
    }

    // Lift cigar
    if (is_first_seg)
        this->lift_cigar(source_contig, aln, start_intvl_idx, pend_start_intvl_idx, num_sclip_start);

    /* Update POS and MPOS */
    if (is_first_seg) {
        if (strand)
            c->pos = c->pos + offset;
        else {
            c->pos = -pos_end + offset + current_intvl.source_start + current_intvl.source_end;
            c->flag ^= BAM_FREVERSE;
        }
    } else {
        if (strand)
            c->mpos = c->mpos + offset;
        else {
            c->mpos = -mpos_end + offset + current_intvl.source_start + current_intvl.source_end;
            c->flag ^= BAM_FMREVERSE;
        }
    }
    // hts_pos_t bam_cigar2qlen(int n_cigar, const uint32_t *cigar);

    return is_liftable;
}


// Pass `dest_contig` by reference because we need it when updating the MD string.
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
        auto mdest_contig = dest_contig;
        bool is_r2_liftable = this->lift_segment(
            aln, hdr, false, mdest_contig, mate_start_intvl_idx, mate_end_intvl_idx);

        // R1 unmapped
        if (c->flag & BAM_FUNMAP) {
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UM_R2_UM;
                // Neither lifted/mapped - set to unmapped
                lo = "UM_UM";
                c->pos = -1;
                c->tid = -1;
                c->qual = 0;
                c->mpos = -1;
                c->mtid = -1;
            } else if (is_r2_liftable) {
                lift_status = LIFT_R1_UM_R2_L;
                // TODO: copy R2 to R1
                c->pos = c->mpos;
                c->tid = c->mtid;
                lo = "UM_L";
            } else {
                lift_status = LIFT_R1_UM_R2_UL;
                // Neither lifted - do nothing
                // c->qual = 255;
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
            c->qual = 255;
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UL_R2_UM;
                // Neither lifted - do nothing
                lo = "UL_UM";
            } else if (is_r2_liftable) {
                lift_status = LIFT_R1_UL_R2_L;
                // #TODO Copy R2 to R1 ??
                c->pos = c->mpos;
                c->tid = c->mtid;
                // c->qual = 255;
                lo = "UL_L";
            } else {
                lift_status = LIFT_R1_UL_R2_UL;
                // Neither lifted - do nothing
                // c->qual = 255;
                lo = "UL_UL";
            }
        }
        // TODO consider STRAND
        c->isize = 
            (lift_status != LIFT_R1_L_R2_L)? 0 : // If any is unliftable/unmapped, set isize to 0
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

ChainMap::ChainMap(std::ifstream& in, int verbose) {
    this->verbose = verbose;
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

// size_t ChainMap::lift_pos(
//     std::string contig, size_t pos,
//     int start_intvl_idx, int end_intvl_idx) {
//     if (start_intvl_idx == -1)
//         return pos;
//     else {
//         auto intvl = this->interval_map[contig][start_intvl_idx];
//         if (intvl.strand) {
//             return pos + intvl.offset;
//         } else {
//             return -pos + intvl.offset + intvl.source_start + intvl.source_end;
//         }
//     }
// }

// Return liftability (bool) and update start/end-interval indexes.
// bool ChainMap::is_liftable(
//     std::string contig, size_t pos,
//     int &start_intvl_idx, int &end_intvl_idx
// ) {
//     int srank = this->get_start_rank(contig, pos);
//     start_intvl_idx = srank - 1;
//     int erank = this->get_end_rank(contig, pos);
//     end_intvl_idx = erank - 1;
// 
//     if (start_intvl_idx == -1 || end_intvl_idx == -1)
//         return false;
//     if (end_intvl_idx - start_intvl_idx == 1)
//         return true;
//     return false;
// }
// 
// std::string ChainMap::lift_contig(
//     std::string contig, int start_intvl_idx, int end_intvl_idx) {
//     if (start_intvl_idx == -1)
//         return "*";
//     else
//         return this->interval_map[contig][start_intvl_idx].target;
// }
//
// void ChainMap::show_interval_info(std::string contig, int pos) {
//     int rank = this->get_start_rank(contig, pos);
//     int intvl_idx = rank - 1;
//     this->interval_map[contig][intvl_idx].debug_print_interval();
// }


size_t ChainMap::lift_pos(
    std::string contig, size_t pos) {
    int rank = this->get_start_rank(contig, pos);
    int intvl_idx = rank - 1;
    if (intvl_idx == -1)
        return pos;
    else {
        auto intvl = this->interval_map[contig][intvl_idx];
        if (intvl.strand) {
            return pos + intvl.offset;
        } else {
            return -pos + intvl.offset + intvl.source_start + intvl.source_end;
        }
    }
}

std::string ChainMap::lift_contig(std::string contig, size_t pos) {
    int rank = this->get_start_rank(contig, pos);
    int intvl_idx = rank - 1;
    if (intvl_idx == -1)
        return "*";
    else
        return this->interval_map[contig][intvl_idx].target;
}

void debug_print_cigar(bam1_t* aln) {
    uint32_t* cigar = bam_get_cigar(aln);
    for (int i = 0; i < aln->core.n_cigar; i++){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        std::cerr << cigar_op_len << bam_cigar_opchr(cigar_op);
    }
    std::cerr << "\n";
}

}
