/*
 * chain.cpp
 *
 * Building a lift-over map from a chain file and performing lift-over
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#include <cmath>
#include <iostream>
#include <queue>
#include "chain.hpp"
#include "leviosam.hpp"
#include "leviosam_utils.hpp"

namespace chain {
Interval::Interval() {
    target = "";
    offset = 0;
    source_start = 0;
    source_end = 0;
    strand = true;
}

Interval::Interval(
    std::string t, int32_t so, int32_t se, int32_t o, bool ss
) {
    target = t;
    offset = o;
    source_start = so;
    source_end = se;
    strand = ss;
}

Interval::Interval(std::ifstream& in) {
    load(in);
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
    size_t str_size;
    in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
    std::vector<char> buf(str_size);
    in.read(reinterpret_cast<char*>(buf.data()), str_size);
    std::string load_target(buf.begin(), buf.end());
    int32_t load_offset;
    in.read(reinterpret_cast<char*>(&load_offset), sizeof(load_offset));
    int32_t load_source_start;
    in.read(reinterpret_cast<char*>(&load_source_start), sizeof(load_source_start));
    int32_t load_source_end;
    in.read(reinterpret_cast<char*>(&load_source_end), sizeof(load_source_end));
    bool load_strand;
    in.read(reinterpret_cast<char*>(&load_strand), sizeof(load_strand));
    target = load_target;
    offset = load_offset;
    source_start = load_source_start;
    source_end = load_source_end;
    strand = load_strand;
}


ChainMap::ChainMap(std::ifstream& in, int verbose, int allowed_intvl_gaps) :
    verbose(verbose), allowed_intvl_gaps(allowed_intvl_gaps){
    load(in);
}


ChainMap::ChainMap(
    std::string fname, int verbose, int allowed_intvl_gaps,
    LengthMap& lm
) : verbose(verbose), allowed_intvl_gaps(allowed_intvl_gaps),
    length_map(lm) {
    std::ifstream chain_f(fname);
    std::string line;
    BitVectorMap start_bv_map;
    BitVectorMap end_bv_map;
    if (chain_f.is_open()) {
        std::string source;
        std::string target;
        int32_t source_offset = 0;
        int32_t target_offset = 0;
        int source_len = 0;
        bool current_ss = true;
        while (getline(chain_f, line)) {
            parse_chain_line(
                line, source, target, source_len, source_offset, 
                target_offset, current_ss,
                start_bv_map, end_bv_map);
        }
        chain_f.close();
    }
    // biv_vector to sd_vector
    for (auto& it: start_bv_map) {
        auto itr = start_map.find(it.first);
        if (itr == start_map.end()) {
            start_map[it.first] = sdsl::sd_vector<>(it.second);
        }
    }
    for (auto& it: end_bv_map) {
        auto itr = end_map.find(it.first);
        if (itr == end_map.end()) {
            end_map[it.first] = sdsl::sd_vector<>(it.second);
        }
    }
    init_rs();
    if (verbose > 2) {
        debug_print_interval_map();
        std::cerr << "\n\nsort intervals\n\n";
    }
    sort_interval_map();
    if (verbose > 2) {
        debug_print_interval_map();
    }
    // TEMP
    interval_map_sanity_check();
}


/* Create start and end bitvectors when see a new `source`.
 * Do nothing if `source` has been seen.
 */
void ChainMap::init_bitvectors(
    std::string source, int source_len,
    BitVectorMap& start_bv_map, BitVectorMap& end_bv_map
) {
    auto itr = start_bv_map.find(source);
    if (itr == start_bv_map.end()) {
        sdsl::bit_vector new_start_bv(source_len), new_end_bv(source_len);
        start_bv_map[source] = new_start_bv;
        end_bv_map[source] = new_end_bv;
    }
}


void ChainMap::sort_interval_map() {
    for (auto &itr: interval_map) {
        sort_intervals(itr.first);
    }
}


/* Sort the intervals in a ChainMap in ascending order */
void ChainMap::sort_intervals(std::string contig) {
    std::sort(
        interval_map[contig].begin(), interval_map[contig].end(),
        [](const Interval& lhs, const Interval& rhs) {
            return lhs.source_start < rhs.source_start;}
    );
}


void ChainMap::debug_print_interval_map() {
    for (auto &intvl_map: interval_map) {
        debug_print_intervals(intvl_map.first);
    }
}


/* Print the intervals in a ChainMap */
void ChainMap::debug_print_intervals(std::string contig) {
    for (auto &intvl: interval_map[contig]) {
        intvl.debug_print_interval();
    }
    std::cerr << "\n";
}


/* Check if the interval map contains any overlaps in the source reference
 * Logic: for each interval, its ending position <= next starting position
 * 
 * Return true if pass; false otherwise.
 */
bool ChainMap::interval_map_sanity_check() {
    for (auto& itr : interval_map){
        std::vector<chain::Interval> v = interval_map[itr.first];
        for (int i = 0; i < v.size()-1; i++) {
            if (v[i].source_end > v[i+1].source_start) {
                std::cerr << "Error: " << itr.first << "\n";
                v[i].debug_print_interval();
                v[i+1].debug_print_interval();
                return false;
            }
        }
    }
    std::cerr << "Interval_map sanity check: passed (no overlaps)\n";
    return true;
}


/* Get the rank in the start bitvector at contig[pos]
 * If pos > len(contig): return rank(contig[-1])
 *
 * Returns:
 *   non-negative int: rank
 *   -1: if pos < 0 or contig does not exist
 */
int ChainMap::get_start_rank(std::string contig, int pos) {
    if (pos < 0) return -1;
    SdVectorMap::const_iterator find_start = start_map.find(contig);
    if (find_start == start_map.end())
        return -1;
    if (pos >= start_rs1_map[contig].size())
        pos = start_rs1_map[contig].size() - 1;
    return start_rs1_map[contig](pos);
}


/* Get the rank in the end bitvector at contig[pos]
 * If pos > len(contig): return rank(contig[-1])
 *
 * Returns:
 *   non-negative int: rank
 *   -1: if pos < 0 or contig does not exist
 */
int ChainMap::get_end_rank(std::string contig, int pos) {
    if (pos < 0)
        return -1;
    SdVectorMap::const_iterator find_end = end_map.find(contig);
    if (find_end == end_map.end())
        return -1;
    if (pos >= end_rs1_map[contig].size())
        pos = end_rs1_map[contig].size() - 1;
    return end_rs1_map[contig](pos);
}

/* Init rank support for all start/end bitvectors */
void ChainMap::init_rs() {
    for (auto& it: start_map) {
        std::string contig = it.first;
        auto itr = start_map.find(contig);
        if (itr == start_map.end()) {
            start_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(start_rs1_map[contig], &start_map[contig]);
    }
    for (auto& it: end_map) {
        std::string contig = it.first;
        auto itr = end_map.find(contig);
        if (itr == end_map.end()) {
            end_rs1_map[contig] = sdsl::sd_vector<>::rank_1_type();
        }
        sdsl::util::init_support(end_rs1_map[contig], &end_map[contig]);
    }
}

/* Parse a chain line and update the ChainMap
 */
void ChainMap::parse_chain_line(
    std::string line, std::string& source, std::string& target,
    int32_t& source_len, int32_t& source_offset,
    int32_t& target_offset, bool& strand,
    BitVectorMap& start_bv_map, BitVectorMap& end_bv_map
) {
    // Split `line` using space as deliminater.
    std::regex s_re("\\s+"); // space
    std::vector<std::string> vec(
        std::sregex_token_iterator(line.begin(), line.end(), s_re, -1),
        std::sregex_token_iterator());
    
    int32_t target_len;
    try{
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
            source_offset = stoi(vec[5]);
            target_offset = (strand)? stoi(vec[10]) : target_len - stoi(vec[10]);
            // length_map.insert(
            //     std::pair<std::string, int32_t>(target, target_len));
            init_bitvectors(source, source_len, start_bv_map, end_bv_map);
        } else if (vec[0] != "") {
            int32_t s_int_start = source_offset;
            int32_t t_int_start = target_offset;
            int32_t s = stoi(vec[0]);
            int32_t s_int_end = source_offset + s;
            int32_t t_int_end = (strand)? target_offset + s : target_offset - s;
            // Set bitvectors
            if (s_int_start > 0)
                start_bv_map[source][s_int_start-1] = 1;
            else if (s_int_start == 0)
                start_bv_map[source][0] = 1;
            if (s_int_end >= source_len)
                end_bv_map[source][source_len-1] = 1;
            else
                end_bv_map[source][s_int_end] = 1;

            if (verbose > 1)
                fprintf(stderr, "source (%d-%d), target (%d-%d)\n",
                        s_int_start, s_int_end, t_int_start, t_int_end);
            int32_t offset = (strand)? t_int_start - s_int_start :
                                       t_int_end - s_int_start;
            Interval intvl(target, s_int_start, s_int_end, offset, strand);
            if (verbose > 1)
                intvl.debug_print_interval();

            IntervalMap::const_iterator find_start =
                interval_map.find(source);
            if (find_start == interval_map.end()) {
                std::vector<chain::Interval> new_intervals;
                interval_map[source] = new_intervals;
            }
            interval_map[source].push_back(intvl);

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
    } catch(...) {
        std::cerr << "Error during parsing chain:\n" << line << "\n";
        exit(1);
    }

    if (verbose > 1) {
        for (auto&& s: vec)
            std::cerr << s << " ";
        std::cerr << "\n";
    }
}


/* Lift CIGAR
 * This version of `lift_cigar()` can be called when no interval indexes are
 * available, but is slower because of additional calls to query for interval indexes.
 */
int ChainMap::lift_cigar(const std::string& contig, bam1_t* aln) {
    bam1_core_t* c = &(aln->core);
    auto pos_end = c->pos + bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));

    // Note that intvl_idx = rank - 1
    auto start_sidx = get_start_rank(contig, c->pos) - 1;
    auto start_eidx = get_end_rank(contig, c->pos) - 1;
    auto end_sidx = get_start_rank(contig, pos_end) - 1;
    auto end_eidx = get_end_rank(contig, pos_end) - 1;

    auto next_intvl = interval_map[contig][start_sidx+1];
    int32_t num_sclip_start = get_num_clipped(
        c->pos, true, contig, start_sidx, start_eidx);
    int32_t num_sclip_end = get_num_clipped(
        pos_end, false, contig, end_sidx, end_eidx);

    return lift_cigar(
        contig, aln, start_sidx, end_sidx,
        num_sclip_start, num_sclip_end);
}


/* Lift CIGAR
 * A fast version of `lift_cigar()`, which is used when auxiliary information 
 * including query start/end interval ranks and #clipped bases is available.
 *
 * The CIGAR string is updated in two cases:
 *   (1) the start-/end-ing position of the segment is outside of an interval, where we set
 *       the overhangs as soft clipped
 *   (2) the segment overlaps one or more gaps between two valid intervals, where we update
 *       the CIGAR with INS or DEL operators
 */
int ChainMap::lift_cigar(
    const std::string& contig, bam1_t* aln,
    const int& start_sidx, const int& end_sidx,
    int num_sclip_start, int num_sclip_end
) {
    if (aln->core.l_qseq > 0 &&
        aln->core.l_qseq < (num_sclip_start + num_sclip_end)){
        return -1;
    }
    // If the read sits perfectly inside an interval, just copy the cigar values
    // (or reverse if in an inversed region):
    //  - no sidx diff
    //  - neither sclip values are positive
    if ((num_sclip_start <= 0) && 
        (num_sclip_end <= 0) &&
        (start_sidx == end_sidx)) {
        uint32_t* cigar = bam_get_cigar(aln);
        if (!interval_map[contig][start_sidx].strand) {
            std::vector<uint32_t> new_cigar;
            for (int i = 0; i < aln->core.n_cigar; i++)
                new_cigar.push_back(cigar[i]);
            for (int i = 0; i < aln->core.n_cigar; i++){
                *(cigar + i) = new_cigar[aln->core.n_cigar - i - 1];
            }
        }
        return 0;
    }

    auto new_cigar = lift_cigar_core(
        contig, aln, start_sidx, end_sidx,
        num_sclip_start, num_sclip_end);

    // There must be >=1 M/I/D/N OPs in the CIGAR, otherwise it's not valid
    bool contain_midn = false;
    for (auto& cg: new_cigar) {
        if (bam_cigar_op(cg) == BAM_CMATCH ||
            bam_cigar_op(cg) == BAM_CINS ||
            bam_cigar_op(cg) == BAM_CDEL ||
            bam_cigar_op(cg) == BAM_CREF_SKIP) {
            contain_midn = true;
            break;
        }
    }
    if (!contain_midn)
        return -1;

    LevioSamUtils::update_cigar(aln, new_cigar);
    return 0;
}


/* Lift one CIGAR run
 */
void ChainMap::lift_cigar_core_one_run(
    std::vector<uint32_t>& new_cigar,
    std::queue<std::tuple<int32_t, int32_t>>& break_points,
    uint32_t cigar_op_len,
    unsigned int cigar_op,
    const uint32_t qlen,
    int& tmp_gap,
    int& query_offset
) {
    int second_half_len = 0;
    // If CIGAR op doesn't consume QUERY, just copy it to `new_cigar`
    if (!(bam_cigar_type(cigar_op) & 1)) {
        push_cigar(new_cigar, cigar_op_len, cigar_op, false);
    } else {
        // Resolve tmp gap bases when we encounter an "M" op
        if (tmp_gap > 0 && cigar_op == BAM_CMATCH) {
            if (tmp_gap < cigar_op_len) {
                cigar_op_len -= tmp_gap;
                tmp_gap = 0;
            } else {
                tmp_gap -= cigar_op_len;
                return;
            }
        }
        auto next_bp = std::get<0>(break_points.front());
        auto next_q_offset = query_offset + cigar_op_len;
        if (verbose >= VERBOSE_DEV) {
            if (tmp_gap != 0)
                std::cerr << "  tmp_gap = " << tmp_gap << "\n";
            std::cerr << "  next_bp = " << next_bp << ", next_q_offset = " << next_q_offset << "\n";
            std::cerr << "  original CIGAR = " << cigar_op_len << bam_cigar_opchr(cigar_op) << "\n";
            if (next_q_offset <= next_bp)
                std::cerr << "  have not reach next_bp: " << cigar_op_len << bam_cigar_opchr(cigar_op) << "\n";
        }
        // Push the current CIGAR OP and advance if have not reached the next breakpoint or
        // there are no remaining breakpoints
        if (next_q_offset <= next_bp || break_points.size() == 0){
            push_cigar(new_cigar, cigar_op_len, cigar_op, false);
            query_offset += cigar_op_len;
        // Split one CIGAR chunk into two parts and insert lift-over bases there
        } else {
            second_half_len = cigar_op_len;
            while (next_q_offset > next_bp && next_bp >= query_offset) {
                auto first_half_len = next_bp - query_offset;
                if (first_half_len > 0) {
                    push_cigar(new_cigar, first_half_len, cigar_op, false);
                    second_half_len -= first_half_len;
                    query_offset += first_half_len;
                }
                int32_t diff = std::get<1>(break_points.front());
                if (verbose >= VERBOSE_DEV) {
                    std::cerr << "  first half: " << first_half_len << bam_cigar_opchr(cigar_op) << " (query_offset -> " << query_offset << ")\n";
                    if (diff > 0)
                        std::cerr << "  update liftover CDEL: " << diff << bam_cigar_opchr(BAM_CDEL) << "\n";
                    else if (diff < 0)
                        std::cerr << "  update liftover CINS: " << -diff << bam_cigar_opchr(BAM_CINS) << "\n";
                }
                // `diff` could be zero when there's a postive-sized poorly 
                // aligned chunk but that will not affect CIGAR updates
                //
                // Unmatched bases in the dest reference
                // (Dest has an insertion w.r.t. source)
                // E.g.
                // read:   TTTT----CCCCCCCGGGG
                // source: TTTT----CCCCCCCGGGG
                // dest:   TTTTAAAACCCCCCCGGGG
                // bp = 4, diff = 4
                // 15M -> 4M4D11M
                if (diff > 0) { // D
                    push_cigar(new_cigar, diff, BAM_CDEL, false);
                // Unmatched bases in the source reference
                // (Dest has an deletion w.r.t source)
                // E.g.
                // read:   TTTTAAAACCCCCCC
                // source: TTTTAAAACCCCCCC
                // dest:   TTTT----CCCCCCC
                // bp = 4, diff = -4
                // 15M -> 4M4I7M
                } else if (diff < 0) { // I
                    // If updating `diff` bases makes the CIGAR longer than QUERY,
                    // truncate `diff` to fit the sequence length.
                    if (query_offset - diff > qlen) {
                        diff = query_offset - qlen;
                    }
                    push_cigar(new_cigar, -diff, BAM_CINS, false);
                    if (cigar_op != BAM_CINS)
                        second_half_len += diff;
                    else
                        tmp_gap -= diff;
                    query_offset -= diff;
                }

                // Popping an empty queue results in undefined behaviors
                if (break_points.size() == 0){
                    break;
                } else {
                    break_points.pop();
                    next_bp = std::get<0>(break_points.front());
                }
            }
            if (second_half_len < 0) {
                // tmp_gap = -second_half_len;
                tmp_gap -= second_half_len;
            } else if (second_half_len > 0) {
                push_cigar(new_cigar, second_half_len, cigar_op, false);
                query_offset += second_half_len;
            }
        }
    }

    if (verbose >= VERBOSE_DEBUG) {
        if (tmp_gap != 0)
            std::cerr << "  tmp_gap = " << tmp_gap << "\n";
        if (second_half_len > 0) {
            std::cerr << "  second_half: " << second_half_len << bam_cigar_opchr(cigar_op) << "\n";
        }
        std::cerr << "  query_offset -> " << query_offset << "\n";
    }
}


/* Core CIGAR lifting function
 *
 * Return a pointer of uint32_t in the htslib CIGAR format.
 */
std::vector<uint32_t> ChainMap::lift_cigar_core(
    const std::string& contig, bam1_t* aln,
    const int& start_sidx, const int& end_sidx,
    const int& num_sclip_start, const int& num_sclip_end
) {
    uint32_t* cigar = bam_get_cigar(aln);
    bam1_core_t* c = &(aln->core);

    // c->l_qseq is someimes zero, e.g secondary alignments
    // We can calculate the actual value by parsing its CIGAR
    uint32_t qlen = (c->l_qseq > 0)?
        c->l_qseq : bam_cigar2qlen(aln->core.n_cigar, cigar);
    if (qlen == 0) {
        std::vector<uint32_t> new_cigar(cigar, cigar + c->n_cigar);
        return new_cigar;
    }

    if (verbose >= VERBOSE_DEBUG) {
        std::cerr << "  num_sclip_start = " << num_sclip_start << "\n";
        std::cerr << "  num_sclip_end = " << num_sclip_end << "\n";
        std::cerr << "    start " << start_sidx << ", pend " << end_sidx << "\n";
        std::cerr << "  * old: ";
        for (int i = 0; i < aln->core.n_cigar; i++){
            auto cigar_op_len = bam_cigar_oplen(cigar[i]);
            auto cigar_op = bam_cigar_op(cigar[i]);
            std::cerr << cigar_op_len << bam_cigar_opchr(cigar_op);
        }
        std::cerr << "\n  pos: " << c->pos << "\n";
    }

    auto break_points = get_bp(contig, c, start_sidx, end_sidx);

    std::vector<uint32_t> new_cigar;
    // Number of bases that need to be clipped in next iterations
    int query_offset = 0;
    int idx = 0;
    // The position of a SAM record starts from the first non-clipped base,
    // so we don't update the first "S" run
    if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) {
        idx += 1;
        auto cigar_op_len = bam_cigar_oplen(cigar[0]);
        push_cigar(new_cigar, cigar_op_len, BAM_CSOFT_CLIP, false);
    }

    int tmp_gap = 0;
    // We first handle soft clipped bases
    if (num_sclip_start > 0) {
        sclip_cigar_front(
            cigar, aln->core.n_cigar, num_sclip_start,
            new_cigar, idx, query_offset);
    }
    for (auto i = idx; i < aln->core.n_cigar; i++){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        if (verbose >= VERBOSE_DEBUG) {
            std::cerr << "\nOP=" << bam_cigar_opchr(cigar[i]) << ", OP_LEN=" << cigar_op_len << "\n";
        }
        if (start_sidx == end_sidx) {
            // If within one interval, update CIGAR and jump to the next CIGAR operator
            push_cigar(new_cigar, cigar_op_len, cigar_op, false);
        } else {
            // Lift one cigar run
            lift_cigar_core_one_run(
                new_cigar, break_points, cigar_op_len, cigar_op,
                qlen, tmp_gap, query_offset);
        }
    }

    // If there are "borrowed" bases (added during liftover) not yet
    // returned, we need to truncate them from the end to make sure
    // CIGAR is valid (bam_cigar2qlen == l_qseq)
    if (tmp_gap > 0) {
        pop_cigar(new_cigar, tmp_gap);
    }

    // Clip from back
    if (num_sclip_end > 0) {
        sclip_cigar_back(new_cigar, num_sclip_end);
    }

    // If there's an "I" in the front or back, replace it with "S"
    if (bam_cigar_op(new_cigar[0]) == BAM_CINS) {
        auto l = bam_cigar_oplen(new_cigar[0]);
        new_cigar[0] = bam_cigar_gen(l, BAM_CSOFT_CLIP);
    }
    if (bam_cigar_op(new_cigar.back()) == BAM_CINS) {
        auto l = bam_cigar_oplen(new_cigar.back());
        new_cigar[new_cigar.size()-1] = bam_cigar_gen(
            l, BAM_CSOFT_CLIP);
    }

    // Reverse CIGAR
    if (!interval_map[contig][start_sidx].strand) {
        std::reverse(new_cigar.begin(), new_cigar.end());
    }

    if (verbose >= VERBOSE_DEBUG) {
        std::cerr << "  len(new_cigar) = " << new_cigar.size() << "\n";
        LevioSamUtils::debug_print_cigar(&new_cigar[0], new_cigar.size());
    }
    return new_cigar;
}


/* Query start_map and end_map for a contig-pos pair and update the indexes
 * Return a bool value:
 *   false if any of the queried contig-pos pair is unavailable
 *   true otherwise
 *
 * sidx: the interval index in start_bv_map[contig]
 * eidx: the interval index in end_bv_map[contig]
 *
 * Note that "interval index = rank - 1"
 */
bool ChainMap::update_interval_indexes(
    const std::string contig, const int32_t pos,
    int32_t& sidx, int32_t& eidx) {
    if (pos < 0) {
        sidx = -1;
        eidx = -1;
        return false;
    }
    // If `contig` is not in the map, set the idx to -1
    // If `pos` is out of scope, query the largest possible position
    SdVectorMap::const_iterator find_start = start_map.find(contig);
    if (find_start == start_map.end()) {
        sidx = -1;
    // TODO
    } else if (pos >= start_rs1_map[contig].size()) {
        auto tmp_pos = start_rs1_map[contig].size() - 1;
        sidx = start_rs1_map[contig](tmp_pos) - 1;
    } else
        sidx = start_rs1_map[contig](pos) - 1;

    SdVectorMap::const_iterator find_end = end_map.find(contig);
    if (find_end == end_map.end()) {
        eidx = -1;
    } else if (pos >= end_rs1_map[contig].size()) {
        auto tmp_pos = end_rs1_map[contig].size() - 1;
        eidx = end_rs1_map[contig](tmp_pos) - 1;
    } else {
        // The first interval: `sidx == 1 && eidx == 0`
        if (end_rs1_map[contig](pos) == 0)
            eidx = 0;
        else
            eidx = end_rs1_map[contig](pos) - 1;
    }


    if (sidx == -1 || eidx == -1)
        return false;
    return true;
}


/* Calculate the distance between a position and nearby chain intervals
 *
 * Return:
 *  - -1 if the query is invalid
 *  - 0 if the position is within a chain interval
 *  - a positive number if there's a positive distance between the position and nearby intervals
 */
int32_t ChainMap::get_num_clipped(
    const int32_t pos, const bool leftmost,
    const std::string& contig, int32_t& sidx, const int32_t& eidx
) {
    if ((sidx <= -1) || (eidx <= -1) || (sidx < eidx)) {
        return -1;
    // Special case: `sidx == 0 && eidx == 0` -> first interval for a contig
    // We can treat this as the `sidx>eidx` case
    } else if (sidx == eidx && eidx != 0) {
        int32_t num_clipped = 0;
        // Advance sidx if we are checking the starting pos of a query
        // Keep sidx unchanged if we are checking the ending pos
        if (leftmost) {
            // If pos is at the left of the rightmost interval, it's not valid
            if (sidx >= interval_map[contig].size() - 1) {
                return -1;
            }
            std::string t = interval_map[contig][sidx].target;
            auto next_intvl = interval_map[contig][sidx+1];
            num_clipped = std::abs(pos - next_intvl.source_start);
            // For a reversed interval, source_start is open
            if (!next_intvl.strand)
                num_clipped += 1;
            sidx += 1;
            // Check if target references are the same
            if (next_intvl.target != t) {
                if (verbose >= VERBOSE_INFO)
                    std::cerr << "[I::chain::get_num_clipped] Mismatched target contig: " <<
                                 t << " and " << interval_map[contig][sidx].target << "\n";
                return -1;
            }
        } else {
            auto curr_intvl = interval_map[contig][sidx];
            num_clipped = std::abs(pos - curr_intvl.source_end);
            // For a forward interval, source_end is open
            if (curr_intvl.strand)
                num_clipped += 1;
        }
        return num_clipped;
    } else if (sidx > eidx || (sidx == eidx && eidx == 0)) {
        std::string t = interval_map[contig][eidx].target;
        for (int i=eidx+1; i <= sidx; i++) {
            if (interval_map[contig][i].target != t) {
                if (verbose >= VERBOSE_INFO)
                    std::cerr << "[I::chain::get_num_clipped] Mismatched target contig: " <<
                                 interval_map[contig][sidx].target <<
                                 " and " << interval_map[contig][eidx].target << "\n";
                return -1;
            }
        }
        return 0;
    }
    std::cerr << "[W::chain::get_num_clipped] Set to unmapped due to an unexpected behavior.\n";
    return -1;
}


/* Lift an alignment segment (there are two segments in a paired-end alignment)
 * Return false if unliftable.
 */
bool ChainMap::lift_segment(
    bam1_t* aln, sam_hdr_t* hdr_source, sam_hdr_t* hdr_dest,
    bool first_seg, std::string& dest_contig
) {
    bam1_core_t* c = &(aln->core);
    // If unmapped, the segment is not liftable.
    // If the aligned REF is not in the map, set the segment to be unliftable.
    if (first_seg) {
        if ((c->flag & BAM_FUNMAP) || (c->tid < 0))
            return false;
    } else {
        if ((c->flag & BAM_FMUNMAP) || (c->mtid < 0))
            return false;
    }
    if (verbose >= VERBOSE_DEBUG) {
        std::cerr << "\n[D::lift_segment] " << bam_get_qname(aln) << "\n";
    }

    int start_sidx = 0;
    int start_eidx = 0;
    auto pos = (first_seg)? c->pos : c->mpos;
    std::string source_contig = (first_seg)?
        hdr_source->target_name[c->tid] : hdr_source->target_name[c->mtid];
    // Update start intervals; if any of the indexes is not valid,
    // the segment is unliftable
    if (!update_interval_indexes(source_contig, pos, start_sidx, start_eidx))
        return false;
    // Set a segment as unliftable if the gap between it
    // and the nearby interval is too large
    //
    // sidx might be advanced here
    int32_t num_sclip_start = get_num_clipped(
        pos, true, source_contig, start_sidx, start_eidx);
    if (verbose >= VERBOSE_DEBUG) {
        // std::cerr << "\n" << bam_get_qname(aln) << "\n";
        debug_print_interval_queries(
            first_seg, true, source_contig, pos, start_sidx, start_eidx);
    }
    if (num_sclip_start < 0) {
        if (verbose >= VERBOSE_INFO) {
            std::cerr << "[I::chain::lift_segment] Set " << bam_get_qname(aln) << 
                         " to unmapped -- invalid num_sclip_start\n";
        }
        return false;
    } else if (num_sclip_start > allowed_intvl_gaps) {
        if (verbose >= VERBOSE_INFO) {
            std::cerr << "[I::chain::lift_segment] Set " << bam_get_qname(aln) << 
                         " to unmapped -- num_sclip_start > " << allowed_intvl_gaps << "\n";
        }
        return false;
    }

    // Lift contig
    auto start_sintvl = interval_map[source_contig][start_sidx];
    dest_contig = ChainMap::lift_contig(start_sintvl);
    if (first_seg)
        c->tid = sam_hdr_name2tid(hdr_dest, dest_contig.c_str());
    else
        c->mtid = sam_hdr_name2tid(hdr_dest, dest_contig.c_str());

    // Check if CONTIG is reasonable
    // Only check mate if paired-end
    if (c->tid < 0 || (c->flag & 1 && c->mtid < 0))
        return false;
    
    // Estimate ending pos of the mate
    // If it's the first segment, we can calculate the query length wrt the REF;
    // If it's the second segment, we can only do a rough estimate b/c we don't know the RLEN of
    // the mate.
    auto pos_end = (first_seg)?
        c->pos + bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln)) :
        c->mpos + c->l_qseq;
    int end_sidx = 0;
    int end_eidx = 0;
    // Update end indexes; if either is unavailable, mark the segment as unliftable
    if (!update_interval_indexes(source_contig, pos_end, end_sidx, end_eidx)){
        return false;
    }
    if (verbose >= VERBOSE_DEBUG) {
        debug_print_interval_queries(
            first_seg, false, source_contig, pos_end, end_sidx, end_eidx);
    }

    int num_sclip_end = 0;
    if (first_seg) {
        num_sclip_end = get_num_clipped(
            pos_end, false, source_contig, end_sidx, end_eidx);
        int num_sclip = num_sclip_start + num_sclip_end;
        if (num_sclip_end < 0) {
            if (verbose >= VERBOSE_INFO) {
                std::cerr << "[I::chain::lift_segment] Set " << bam_get_qname(aln) << 
                             " to unmapped -- invalid num_sclip_end\n";
            }
            return false;
        }
        // #SCLIP must <= `allowed_intvl_gaps`
        if (num_sclip > allowed_intvl_gaps) {
            if (verbose >= VERBOSE_INFO) {
                std::cerr << "[I::chain::lift_segment] Set " << bam_get_qname(aln) << 
                             " to unmapped -- (num_sclip_start+num_sclip_end) > " <<
                             allowed_intvl_gaps << "\n";
            }
            return false;
        }
    }

    // If the corresponding intervals of the left and right ends differ,
    // check if the intervals are legal.
    // Intervals with too large gap, inconsistent target/strandness will
    // not pass the check.
    if (start_sidx != end_sidx) {
        // Returns true if legal
        if (!check_multi_intvl_legality(
            source_contig, aln, start_sidx, end_sidx, allowed_intvl_gaps)
        ) {
            return false;
        }
    }

    auto next_sintvl = interval_map[source_contig][end_sidx];
    // Lift CIGAR
    auto rlen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));
    if (first_seg && lift_cigar(
        source_contig, aln, start_sidx, end_sidx,
        num_sclip_start, num_sclip_end) != 0) {
        if (verbose >= VERBOSE_INFO) {
            std::cerr << "[I::chain::lift_segment] Set " << bam_get_qname(aln) << 
                         " to unmapped b/c of failure in lift_cigar\n";
        }
        return false;
    }

    if (!(c->flag & BAM_FUNMAP) && !(c->flag & BAM_FSECONDARY) && !(c->flag & BAM_FSUPPLEMENTARY)) {
        uint32_t* cigar = bam_get_cigar(aln);
        if (c->l_qseq != bam_cigar2qlen(aln->core.n_cigar, cigar)) {
            std::cerr << "[W::chain::lift_segment] Read " << bam_get_qname(aln) << " (flag=" << 
                         c->flag << " has mismatched c->l_qseq (" << c->l_qseq << 
                         ") and bam_cigar2qlen (" << bam_cigar2qlen(aln->core.n_cigar, cigar) << 
                         "). This is not an expected behavior\n";
            LevioSamUtils::debug_print_cigar(cigar, c->n_cigar);
            return false;
        }
    }

    // Lift POS and MPOS
    // Updating positions affects the lift-over of other fields (e.g. CIGAR), so this has to be
    // performed in the end of the lift-over process.
    if (start_sintvl.strand) {
        lift_pos(aln, pos_end, start_sintvl, first_seg);
        if (first_seg && num_sclip_start > 0)
            c->pos += num_sclip_start;
    } else {
        lift_pos(aln, pos_end, next_sintvl, first_seg);
        if (first_seg && num_sclip_end > 0)
            c->pos += num_sclip_end;
    }

    // There can be cases where a position is lifted to a negative value,
    // in which we set the read to unmapped.
    if (c->pos < 0) {
        std::cerr << "[W::chain::lift_segment] Read " << bam_get_qname(aln) << 
                     " is lifted to a negative position (" << c->pos << 
                     "). Set it to unmapped\n";
        if (verbose >= VERBOSE_DEBUG) {
            std::cerr << "flag=" << c->flag << "; num_sclip_start=" << num_sclip_start <<
                         ", num_sclip_end=" << num_sclip_end << "\n";
            LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), c->n_cigar);
        }
        c->pos = 0;
        return false;
    } else if ((c->flag & BAM_FPAIRED) && c->mpos < 0) {
        // Invalid mate positions can be benign since there's no CIGAR information to add
        // SOFT_CLIP based in the beginning. We report these as INFO but not WARNING
        if (verbose >= VERBOSE_INFO) {
            std::cerr << "[I::chain::lift_segment] The mate of read " << bam_get_qname(aln) << 
                         " is lifted to a negative position (" << c->mpos << ")\n";
        }
        c->mpos = 0;
        return false;
    }
    return true;
}


// Pass `dest_contig` by reference because we need it when updating the MD string.
void ChainMap::lift_aln(
    bam1_t* aln, sam_hdr_t* hdr_source, sam_hdr_t* hdr_dest,
    std::string& dest_contig
) {
    bam1_core_t* c = &(aln->core);
    uint16_t flag = c->flag;
    hts_pos_t pos = c->pos;
    hts_pos_t mpos = c->mpos;

    bool r1_liftable = lift_segment(aln, hdr_source, hdr_dest, true, dest_contig);
    if (!r1_liftable) {
        LevioSamUtils::update_flag_unmap(aln, true);
    }

    size_t lift_status;
    std::string lo;

    // Paired
    if (c->flag & BAM_FPAIRED) {
        std::string null_mdest_contig;
        bool r2_liftable = lift_segment(aln, hdr_source, hdr_dest, false, null_mdest_contig);
        if (!r2_liftable)
            LevioSamUtils::update_flag_unmap(aln, false);

        // We currently don't use the "unliftable" category - 
        // an unliftable read is considered as unmapped.
        // R1 unmapped
        if (c->flag & BAM_FUNMAP) {
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UM_R2_UM;
                // Neither lifted/mapped - set to unmapped
                lo = "UM_UM";
                c->pos = -1;
                c->tid = -1;
                c->mpos = -1;
                c->mtid = -1;
            } else if (r2_liftable) {
                lift_status = LIFT_R1_UM_R2_L;
                // Copy R2 to R1
                c->pos = c->mpos;
                c->tid = c->mtid;
                lo = "UM_L";
            } else {
                lift_status = LIFT_R1_UM_R2_UL;
                // Neither is lifted
                lo = "UM_UL";
            }
        // R1 mapped, liftable
        } else if (r1_liftable) {
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_L_R2_UM;
                // Copy r1 to r2
                c->mpos = c->pos;
                c->mtid = c->tid;
                lo = "L_UM";
            } else if (r2_liftable) {
                lift_status = LIFT_R1_L_R2_L;
                lo = "L_L";
            } else {
                lift_status = LIFT_R1_L_R2_UL;
                // Copy R1 to R2
                c->mpos = c->pos;
                c->mtid = c->tid;
                lo = "L_UL";
            }
        // R1 mapped, un-liftable
        } else {
            // c->qual = 255;
            if (c->flag & BAM_FMUNMAP) {
                lift_status = LIFT_R1_UL_R2_UM;
                // Neither is lifted
                lo = "UL_UM";
            } else if (r2_liftable) {
                lift_status = LIFT_R1_UL_R2_L;
                // Copy R2 to R1
                c->pos = c->mpos;
                c->tid = c->mtid;
                lo = "UL_L";
            } else {
                lift_status = LIFT_R1_UL_R2_UL;
                // Neither is lifted
                lo = "UL_UL";
            }
        }
        // Lift the TLEN (c->isize) field
        if ((lift_status != LIFT_R1_L_R2_L) || // If any is unliftable/unmapped
            (c->isize == 0) || // The reads were originally not paired properly
            (c->tid != c->mtid)) // Mismatched contigs
            // !(c->flag & 2)) // Not properly paired
        { 
            c->isize = 0;
        } else {
            if (flag == c->flag)
                c->isize = c->isize + (mpos - pos) - (c->mpos - c->pos);
            else
                c->isize = -(c->isize + (mpos - pos) + (c->mpos - c->pos));
        }

    // Unpaired
    } else {
        if (c->flag & BAM_FUNMAP) {
            lift_status = LIFT_R_UM;
            lo = "UM";
        } else if (r1_liftable) {
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
    for (auto& x: interval_map) {
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
    nelems = start_map.size();
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    size += sizeof(nelems);
    for (auto& x: start_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        size += sizeof(str_size) + str_size;
        size += x.second.serialize(out);
    }
    // end_map
    std::cerr << "serializing end_maps...\n";
    nelems = end_map.size();
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    size += sizeof(nelems);
    for (auto& x: end_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        size += sizeof(str_size) + str_size;
        size += x.second.serialize(out);
    }

    // length_map
    std::cerr << "serializing length_maps...\n";
    size += LevioSamUtils::serialize_lengthmap(out, length_map);
    return size;
}


// loads from stream
void ChainMap::load(std::ifstream& in) {
    size_t map_size;
    if (verbose > 2)
        std::cerr << "Reading interval map...\n";
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
            interval_map[key].push_back(Interval(in));
        }
    }
    if (verbose > 2)
        std::cerr << "Reading start BV...\n";
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        sdsl::sd_vector<> sdv;
        sdv.load(in);
        start_map[key] = sdv;
    }
    if (verbose > 2)
        std::cerr << "Reading end BV...\n";
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        sdsl::sd_vector<> sdv;
        sdv.load(in);
        end_map[key] = sdv;
    }
    if (verbose > 2)
        std::cerr << "Reading length map...\n";
    length_map = LevioSamUtils::load_lengthmap(in);
    if (verbose > 2)
        std::cerr << "Initialize BVs...\n";
    init_rs();
}


/* Lift over the POS field of an alignment
 *   - aln: target aln object
 *   - pos_end: putated ending position of an aln
 *   - intvl: Interval object of the target
 *   - first_seg: true if is read1/single-end; false if read2
 */
void ChainMap::lift_pos(
    bam1_t* aln, const hts_pos_t& pos_end,
    const chain::Interval& intvl, const bool& first_seg
) {
    auto c = &(aln->core);
    if (first_seg) {
        if (intvl.strand) {
            c->pos = lift_pos(c->pos, intvl);
        } else {
            c->pos = lift_pos(pos_end, intvl);
            // Update FLAG, SEQ and QUAL if in a reversed chain
            c->flag ^= BAM_FREVERSE;
            LevioSamUtils::reverse_seq_and_qual(aln);
        }
    } else {
        if (intvl.strand) {
            c->mpos = ChainMap::lift_pos(c->mpos, intvl);
        } else {
            c->mpos = lift_pos(pos_end, intvl);
            // Update FLAG if in a reversed chain
            c->flag ^= BAM_FMREVERSE;
        }
    }
}


// This saves one `rank` query by providing `intvl_idx`
hts_pos_t ChainMap::lift_pos(
    const hts_pos_t& pos, const chain::Interval& intvl
) {
    if (intvl.strand) {
        return pos + intvl.offset;
    } else {
        return -pos + intvl.offset + intvl.source_start + intvl.source_end;
    }
}


/* Lift over a position
 *
 * Inputs:
 *  - contig
 *  - pos
 *  - allowed_gaps
 *  - left: true if the POS is at the left side of an interval; otherwise false
 */
hts_pos_t ChainMap::lift_pos(
    const std::string& contig, const hts_pos_t& pos,
    const int& allowed_gaps, const bool& left
) {
    int sidx = 0;
    int eidx = 0;
    if (!update_interval_indexes(contig, pos, sidx, eidx)) {
        return -1;
    }
    if (left) {
        int32_t num_sclip_start = get_num_clipped(
            pos, true, contig, sidx, eidx);
        if (num_sclip_start > allowed_gaps) {
            return -1;
        }
        chain::Interval intvl = this->interval_map[contig][sidx];
        return lift_pos(pos, intvl) + num_sclip_start;
    } else {
        int32_t num_sclip_end = get_num_clipped(
            pos, false, contig, sidx, eidx);
        if (num_sclip_end > allowed_gaps) {
            return -1;
        }
        chain::Interval intvl = this->interval_map[contig][sidx];
        return lift_pos(pos, intvl) - num_sclip_end;
    }
    return -1;
}


std::string ChainMap::lift_contig(const chain::Interval& intvl) {
    return intvl.target;
}


std::string ChainMap::lift_contig(
    const std::string& contig, const hts_pos_t& pos) {
    int rank = this->get_start_rank(contig, pos);
    int intvl_idx = rank - 1;
    if (intvl_idx == -1)
        return "*";
    else {
        auto intvl = this->interval_map[contig][intvl_idx];
        return lift_contig(intvl);
    }
}


/* Update a cigar vector
 * It is essentially equivalent with 
 * `cigar.push_back(bam_cigar_gen(len, op));`
 * 
 * The difference is this function additionally performs 
 * CIGAR operator reduction, e.g.
 *   (1M2M -> 3M), (1D1I -> 1M)
 * 
 * We set `no_reduct=true` to disable the reduction functionality.
 */
void push_cigar(
    std::vector<uint32_t>& cigar, uint32_t len,
    uint16_t op, const bool no_reduct=false
) {
    if (len == 0)
        return;
    if (cigar.size() == 0 || no_reduct) {
        cigar.push_back(bam_cigar_gen(len, op));
        return;
    }
    auto back_op = bam_cigar_op(cigar.back());
    auto back_type = bam_cigar_type(back_op);
    auto back_len = bam_cigar_oplen(cigar.back());
    auto op_type = bam_cigar_type(op);
    // If operators are the same, merge
    // We also merge S-I/I-S cases
    // We don't merge N-D, H-P because they might mean differently
    if (back_op == op || (back_type == 1 && op_type == 1)) {
        len += back_len;
        cigar.back() = bam_cigar_gen(len, back_op);
    // Cancel out complementary operators
    } else if ((back_type == 2 && op_type == 1) ||
               (back_type == 1 && op_type == 2)) {
        // For cases S-D, we ignore Ds, since deletions don't 
        // contribute to either the alignment or the position 
        // in soft-clipped regions.
        if (back_op == BAM_CSOFT_CLIP) {
            return;
        } if (len == back_len) {
            cigar.pop_back();
            push_cigar(cigar, len, BAM_CMATCH, false);
        } else if (len > back_len){
            len -= back_len;
            cigar.back() = bam_cigar_gen(back_len, BAM_CMATCH);
            cigar.push_back(bam_cigar_gen(len, op));
        } else {
            back_len -= len;
            cigar.back() = bam_cigar_gen(back_len, back_op);
            cigar.push_back(bam_cigar_gen(len, BAM_CMATCH));
        }
    } else
        cigar.push_back(bam_cigar_gen(len, op));
}


/* Pop `size` bases from a CIGAR (represented as a vector)
 */
void pop_cigar(std::vector<uint32_t>& cigar, uint32_t size) {
    while (size > 0) {
        auto cg = cigar.back();
        cigar.pop_back();
        auto cigar_op_len = bam_cigar_oplen(cg);
        auto cigar_op = bam_cigar_op(cg);
        // Check if the last operator consumes REF; if it doesn't, just pop it
        if (bam_cigar_type(cigar_op) & 1) {
            if (cigar_op_len > size) {
                cigar_op_len -= size;
                size = 0;
                push_cigar(cigar, cigar_op_len, cigar_op, false);
            } else if (cigar_op_len == size) {
                size = 0;
            } else {
                size -= cigar_op_len;
            }
        }
    }
}


void sclip_cigar_front(
    uint32_t* cigar, const uint32_t& n_cigar, int len_clip,
    std::vector<uint32_t>& new_cigar, int& idx, int& query_offset
) {
    if (len_clip == 0)
        return;

    for (auto i = idx; i < n_cigar; i++){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        // Only replace bases that consume the QUERY with the SOFT_CLIP ("S") operator
        if (bam_cigar_type(cigar_op) & 1) {
            // If a CIGAR OP is longer than #clipped, append clipped bases to the updated CIGAR
            // and truncate the original CIGAR for future usage.
            if (cigar_op_len >= len_clip) {
                cigar_op_len -= len_clip;
                push_cigar(new_cigar, len_clip, BAM_CSOFT_CLIP, false);
                query_offset += len_clip;
                // Update `idx` so that we don't need to revisit
                // previous CIGAR OPs.
                if (cigar_op_len > 0) {
                    cigar[i] = bam_cigar_gen(cigar_op_len, cigar_op);
                    idx = i;
                } else {
                    idx = i + 1;
                }
                return;
            } else {
                push_cigar(new_cigar, cigar_op_len, BAM_CSOFT_CLIP, false);
                query_offset += cigar_op_len;
                len_clip -= cigar_op_len;
            }
        }
    }
}


void sclip_cigar_back(
    std::vector<uint32_t>& cigar, int len_clip
) {
    if (len_clip == 0)
        return;

    int remaining_lc = len_clip;
    for (int i = cigar.size() - 1; i >= 0; i--){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        cigar.pop_back();
        // Only replace bases that consume the QUERY with the SOFT_CLIP ("S") operator
        if (bam_cigar_type(cigar_op) & 1) {
            if (cigar_op_len >= remaining_lc) {
                // If a CIGAR OP is longer than #clipped, append clipped bases
                // to the updated CIGAR and truncate the original CIGAR.
                cigar_op_len -= remaining_lc;
                push_cigar(cigar, cigar_op_len, cigar_op, false);
                push_cigar(cigar, len_clip, BAM_CSOFT_CLIP, false);
                return;
            } else if (cigar_op == BAM_CMATCH) {
                // If a CIGAR OP is not long enough, shorten `remaining_lc`
                // and proceed to an earlier OP
                remaining_lc -= cigar_op_len;
            } else if (cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
                len_clip += cigar_op_len;
            } else {
                std::cerr << "[W::chain::sclip_cigar_back] Unexpected OP during back_clipping: " << cigar_op << "\n";
            }
        }
    }
}


std::queue<std::tuple<int32_t, int32_t>> ChainMap::get_bp(
    const std::string& contig, const bam1_core_t* const c,
    const int& start_sidx, const int& end_sidx
){
    std::queue<std::tuple<int32_t, int32_t>> break_points;
    for (auto i = start_sidx; i <= end_sidx; i++) {
        int bp = interval_map[contig][i].source_end - c->pos;
        int diff = interval_map[contig][i+1].offset -
                   interval_map[contig][i].offset;
        if (interval_map[contig][i+1].target != interval_map[contig][i].target)
            continue;
        if (bp > 0) {
            break_points.push(std::make_tuple(bp, diff));
            // DEBUG
            if (verbose >= VERBOSE_DEBUG) {
                interval_map[contig][i].debug_print_interval();
                std::cerr << "bp=" << bp << ", diff=" << diff << "\n";
            }
            // END_DEBUG
        }
    }
    return break_points;
}


/* If an alignment is overlapped with multiple intervals, check the gap size between the 
 * intervals. If either
 *   (1) the offset difference between any adjacent interval pair 
 *   (2) the total offset difference between the first and the last interval 
 * is greater than `allowed_intvl_gaps`, mark the alignment as unliftable (return false).
 */
bool ChainMap::check_multi_intvl_legality(
    const std::string& s, bam1_t* aln,
    const int& start_sidx, const int& end_sidx,
    const int& allowed_intvl_gaps
) {
    Interval start_sintvl = interval_map[s][start_sidx];
    Interval next_sintvl = interval_map[s][end_sidx];
    // Check the gap between the first and the last intervals
    if (std::abs(next_sintvl.offset - start_sintvl.offset) > allowed_intvl_gaps) {
        if (verbose >= VERBOSE_INFO) {
            std::cerr << "[I::chain::check_multi_intvl_legality] Set " << bam_get_qname(aln) << 
                         " to unmapped -- 1st-and-last chain gap > " << allowed_intvl_gaps << "\n";
        }
        return false;
    }
    // Check all interval pairs
    for (auto j = start_sidx; j < end_sidx; j ++) {
        if (interval_map[s][j].strand != interval_map[s][j+1].strand) {
            if (verbose >= VERBOSE_INFO) {
                std::cerr << "[I::chain::check_multi_intvl_legality] Set " << bam_get_qname(aln) << 
                             " to unmapped -- intervals not having concordant strandness\n";
            }
            return false;
        }
        if (interval_map[s][j].target != interval_map[s][j+1].target) {
            if (verbose >= VERBOSE_INFO) {
                std::cerr << "[I::chain::check_multi_intvl_legality] Set " << bam_get_qname(aln) << 
                             " to unmapped -- intervals not having concordant targets\n";
            }
            return false;
        }
        if (std::abs(interval_map[s][j+1].offset -
                     interval_map[s][j].offset) > allowed_intvl_gaps) {
            if (verbose >= VERBOSE_INFO) {
                std::cerr << "[I::chain::check_multi_intvl_legality] Set " << bam_get_qname(aln) << 
                             " to unmapped -- >=1 chain gaps > " << allowed_intvl_gaps << "\n";
            }
            return false;
        }
    }
    return true;
}


void ChainMap::debug_print_interval_queries(
    const bool first_seg, const bool leftmost,
    const std::string contig, const int32_t pos,
    const int32_t sidx, const int32_t eidx) {
    if (first_seg) {
        std::cerr << "First segment ";
        if (leftmost) std::cerr << "(start)\n";
        else std::cerr << "(estimated end)\n";
    } else {
        std::cerr << "Second segment";
        if (leftmost) std::cerr << "(start)\n";
        else std::cerr << "(estimated end)\n";
    }
    std::cerr << "  Source: " << contig << ":" << pos << "\n";
    std::cerr << "  Queried interval indexes:\n";
    std::cerr << "    start: " << sidx << " ";
    interval_map[contig][sidx].debug_print_interval();
    std::cerr << "    end  : " << eidx << " ";
    interval_map[contig][eidx].debug_print_interval();
}

}
