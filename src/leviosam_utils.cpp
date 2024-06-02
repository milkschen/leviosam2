/*
 * leviosam_utils.cpp
 *
 * Utility functions for the leviosam software
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "leviosam_utils.hpp"

#include <algorithm>
#include <cmath>
#include <regex>

size_t kstr_get_m(size_t var) {
    size_t lvar = (size_t)std::exp2(std::ceil(std::log2(var)));
    return lvar;
}

namespace LevioSamUtils {

void WriteDeferred::init(const SplitRules& split_rules, sam_hdr_t* ihdr,
                         sam_hdr_t* ohdr, const BedUtils::Bed& b_defer_source,
                         const BedUtils::Bed& b_defer_dest,
                         const BedUtils::Bed& b_commit_source,
                         const BedUtils::Bed& b_commit_dest,
                         const float& b_isec_th) {
    write_deferred = true;
    for (auto& r : split_rules) {
        split_modes.emplace(r.first);
        if (r.first == "mapq") {
            min_mapq = r.second;
        } else if (r.first == "clipped_frac") {
            max_clipped_frac = r.second;
        } else if (r.first == "isize") {
            max_isize = r.second;
        } else if (r.first == "aln_score") {
            min_aln_score = r.second;
        } else if (r.first == "hdist") {
            max_hdist = r.second;
        }
    }
    hdr = ohdr;
    hdr_orig = ihdr;
    bed_defer_source = b_defer_source;
    bed_defer_dest = b_defer_dest;
    bed_commit_source = b_commit_source;
    bed_commit_dest = b_commit_dest;
    bed_isec_threshold = b_isec_th;
}

/** Closes the deferred and unliftable HTS files. */
void WriteDeferred::close_sams() {
    if (write_deferred) {
        sam_close(out_fp);
        sam_close(out_fp_orig);
    } else {
        std::cerr << "[W::WriteDeferred::close_sams] No HTS files to close "
                     "when `write_deferred` is False.\n";
    }
};

void WriteDeferred::print_info() {
    std::cerr << "[I::WriteDeferred::print_info] Alignments with any of the "
                 "below features are deferred:\n";
    std::cerr << " - unlifted\n";
    if (split_modes.find("mapq") != split_modes.end()) {
        std::cerr << " - MAPQ (pre-liftover) < " << min_mapq << "\n";
    }
    if (split_modes.find("isize") != split_modes.end()) {
        std::cerr << " - TLEN/isize (post-liftover) > " << max_isize << "\n";
        std::cerr << " - TLEN/isize (post-liftover) == 0\n";
    }
    if (split_modes.find("clipped_frac") != split_modes.end()) {
        std::cerr << " - Fraction of clipped bases (post-liftover) >= "
                  << max_clipped_frac << "\n";
    }
    if (split_modes.find("aln_score") != split_modes.end()) {
        std::cerr << " - AS:i (pre-liftover) < " << min_aln_score << "\n";
    }
    if (split_modes.find("hdist") != split_modes.end()) {
        std::cerr << " - NM:i (post-liftover) > " << max_hdist << "\n";
    }
    if (bed_defer_source.get_fn().size() > 0) {
        std::cerr << " - In BED deferred (source) " << bed_defer_source.get_fn()
                  << "\n";
    }
    if (bed_defer_dest.get_fn().size() > 0) {
        std::cerr << " - In BED deferred (dest) " << bed_defer_dest.get_fn()
                  << "\n";
    }
    if (bed_commit_source.get_fn().size() > 0) {
        std::cerr << " - Not in BED commit (source) "
                  << bed_commit_source.get_fn() << "\n";
    }
    if (bed_commit_dest.get_fn().size() > 0) {
        std::cerr << " - Not in BED commit (dest) " << bed_commit_dest.get_fn()
                  << "\n";
    }
}

/** Gets the write_deferred status. */
bool WriteDeferred::get_write_deferred() { return write_deferred; }

/** Creates deferred and unliftable HTS files.
 *
 * @param out_prefix Output prefix
 * @param out_format Output format. Options: "sam", "bam"
 */
void WriteDeferred::open_sams(const std::string out_prefix,
                              const std::string out_format) {
    std::string out_mode = (out_format == "sam") ? "w" : "wb";
    std::string out_fn = out_prefix + "-deferred." + out_format;
    out_fp = sam_open(out_fn.data(), out_mode.data());
    if (sam_hdr_write(out_fp, hdr) < 0) {
        std::cerr << "[E::WriteDeferred::init] Failed to write sam_hdr for "
                  << out_fn << "\n";
        exit(1);
    }

    std::string out_fn_orig = out_prefix + "-unliftable." + out_format;
    out_fp_orig = sam_open(out_fn_orig.data(), out_mode.data());
    if (sam_hdr_write(out_fp_orig, hdr_orig) < 0) {
        std::cerr << "[E::WriteDeferred::init] Failed to write sam_hdr for "
                  << out_fn_orig << "\n";
        exit(1);
    }
}

/* Returns true if an alignment is excluded (committed)
 */
bool WriteDeferred::commit_aln_source(const bam1_t* const aln) {
    // If defer mode is not activated, all reads are committed
    if (!write_deferred) return true;

    const bam1_core_t* c = &(aln->core);
    if (c->flag & BAM_FUNMAP) return false;

    std::string rname = hdr_orig->target_name[c->tid];
    auto rlen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));
    size_t pos_end = c->pos + rlen;
    // The commit regions need to be considered ahead of other rules,
    // i.e. a low MAPQ read in a commit region should be committed
    if (bed_commit_source.intersect(rname, c->pos, pos_end,
                                    bed_isec_threshold)) {
        return true;
    }
    if (bed_defer_source.intersect(rname, c->pos, pos_end,
                                   bed_isec_threshold)) {
        return false;
    }
    return false;
}

/* Returns true if an alignment is committed
 */
bool WriteDeferred::commit_aln_dest(const bam1_t* const aln) {
    // If defer mode is not activated, all reads are committed
    if (!write_deferred) return true;

    const bam1_core_t* c = &(aln->core);

    if (c->flag & BAM_FUNMAP) return false;
    // MAPQ and AS are pre-lift-over values, but we need to put the logic here
    // since our commit-defer decision is based on (commit_aln_source |
    // commit_aln_dest)
    if (split_modes.find("mapq") != split_modes.end()) {
        if (c->qual < min_mapq) return false;
    }
    if (split_modes.find("aln_score") != split_modes.end()) {
        if (bam_aux2i(bam_aux_get(aln, "AS")) < min_aln_score) return false;
    }
    // End MAPQ and AS
    if (split_modes.find("isize") != split_modes.end()) {
        if (c->isize == 0 || c->isize > max_isize || c->isize < -max_isize)
            return false;
    }

    std::string rname = hdr->target_name[c->tid];
    auto rlen = bam_cigar2rlen(c->n_cigar, bam_get_cigar(aln));
    size_t pos_end = c->pos + rlen;
    if (bed_defer_dest.intersect(rname, c->pos, pos_end, bed_isec_threshold)) {
        return false;
    }
    if (bed_commit_dest.intersect(rname, c->pos, pos_end, bed_isec_threshold)) {
        return true;
    }
    if (split_modes.find("clipped_frac") != split_modes.end()) {
        if (1.0 - (rlen / c->l_qseq) > max_clipped_frac) return false;
    }
    if (split_modes.find("hdist") != split_modes.end()) {
        if (bam_aux2i(bam_aux_get(aln, "NM")) > max_hdist) return false;
    }
    return true;
}

/** Removes the MN:i and MD:z tags from an alignment object.
 *
 * @param aln Alignment object.
 */
void remove_nm_md_tag(bam1_t* aln) {
    uint8_t* ptr = NULL;
    if ((ptr = bam_aux_get(aln, "MD")) != NULL) {
        bam_aux_del(aln, ptr);
    }
    if ((ptr = bam_aux_get(aln, "NM")) != NULL) {
        bam_aux_del(aln, bam_aux_get(aln, "NM"));
    }
}

void remove_aux_tag(bam1_t* aln, const std::string tag) {
    uint8_t* ptr = NULL;
    if ((ptr = bam_aux_get(aln, tag.c_str())) != NULL) {
        bam_aux_del(aln, ptr);
    }
}

/** Gets the reference length of the mate segment.
 *
 * Leverages the "MC:Z" flag when available. If not, infers the query length.
 *
 * @param aln Alignment object.
 * @return Reference length of the mate segment.
 */
hts_pos_t get_mate_query_len_on_ref(const bam1_t* aln) {
    const bam1_core_t* c = &(aln->core);
    uint8_t* ptr = NULL;
    uint8_t* mc = bam_aux_get(aln, "MC");
    if (mc) {
        char* mate_cigar = bam_aux2Z(mc);
        uint32_t* a_cigar = NULL;
        size_t n_cigar = 0;
        sam_parse_cigar(mate_cigar, NULL, &a_cigar, &n_cigar);
        return bam_cigar2rlen(n_cigar, a_cigar);
    } else {
        return c->l_qseq;
    }
}

/** Gets the forward read sequence.
 *
 * Adapted from: https://github.com/samtools/samtools/blob/develop/bam_fastq.c
 *
 * @param rec Alignment obejct.
 * @return Read sequence string.
 */
std::string get_forward_read(const bam1_t* rec) {
    char* seq = (char*)bam_get_seq(rec);
    std::string read = "";

    for (int n = 0; n < rec->core.l_qseq; n++) {
        if (rec->core.flag & BAM_FREVERSE)
            read.append(1, seq_nt16_str[seq_comp_table[bam_seqi(seq, n)]]);
        else
            read.append(1, seq_nt16_str[bam_seqi(seq, n)]);
    }
    if (rec->core.flag & BAM_FREVERSE) std::reverse(read.begin(), read.end());
    return read;
}

/** Gets the aligned read sequence.
 *
 * Use `get_forward_read()` if interested in the forward read sequence.
 *
 * @param rec Alignment obejct.
 * @return Read sequence string.
 */
std::string get_read(const bam1_t* rec) {
    char* seq = (char*)bam_get_seq(rec);
    std::string read = "";

    for (int n = 0; n < rec->core.l_qseq; n++) {
        read.append(1, seq_nt16_str[bam_seqi(seq, n)]);
    }
    return read;
}

/**
 * Sets an alignment to unmapped.
 *   - Clear forward/reverse status.
 *   - If paired, changed to improper paired.
 *
 * An unliftable read is considered as unmapped. The BAM_FUNMAP or BAM_FMUNMAP
 * flags are updated accordingly.
 *
 * @param aln
 * @param first_set
 * @param keep_mapq
 */
void update_flag_unmap(bam1_t* aln, const bool first_seg,
                       const bool keep_mapq = false) {
    bam1_core_t* c = &(aln->core);
    if (first_seg) {
        if (c->flag & BAM_FREVERSE) LevioSamUtils::reverse_seq_and_qual(aln);
        c->flag |= BAM_FUNMAP;
        c->flag &= ~BAM_FPROPER_PAIR;
        c->flag &= ~BAM_FREVERSE;
        if (!keep_mapq) c->qual = 0;
        c->pos = -1;
        c->tid = -1;
    } else {
        c->flag |= BAM_FMUNMAP;
        c->flag &= ~BAM_FPROPER_PAIR;
        c->flag &= ~BAM_FMREVERSE;
    }
    // Keep the secondary/supplementary annotation if unmapped
    // This might violate some strict SAM linter, but keeping the
    // annotations can avoid converting SUPP reads to FASTQ
    // c->flag &= ~BAM_FSECONDARY;
    // c->flag &= ~BAM_FSUPPLEMENTARY;
}

void WriteDeferred::write_deferred_bam(bam1_t* aln, sam_hdr_t* hdr) {
    if ((aln->core.flag &
         BAM_FSECONDARY) ||                // Secondary alignment - no SEQ field
        (aln->core.flag & BAM_FQCFAIL) ||  // not passing filters
        (aln->core.flag & BAM_FDUP) ||     // PCR or optinal duplicate
        (aln->core.flag & BAM_FSUPPLEMENTARY)) {  // supplementary alignment
        return;
    }
    auto w_ret = sam_write1(out_fp, hdr, aln);
    if (w_ret < 0) {
        std::cerr
            << "[E::WriteDeferred::write_deferred_bam] Failed to write record "
            << bam_get_qname(aln) << "\n";
        exit(1);
    }
}

void WriteDeferred::write_deferred_bam_orig(bam1_t* aln) {
    auto w_ret = sam_write1(out_fp_orig, hdr_orig, aln);
    if (w_ret < 0) {
        std::cerr << "[E::WriteDeferred::write_deferred_bam_orig] Failed to "
                     "write record "
                  << bam_get_qname(aln) << "\n";
        exit(1);
    }
}

/* Construct a FastqRecord object from seq and qual */
FastqRecord::FastqRecord(const std::string& seq, const std::string& qual) {
    seq_str = seq;
    qual_str = qual;
}

/* Construct a FastqRecord object from a BAM record */
FastqRecord::FastqRecord(bam1_t* a) {
    aln = bam_init1();
    if (bam_copy1(aln, a) == NULL) {
        std::cerr << "Failed to copy " << bam_get_qname(a) << "\n";
    };
}

FastqRecord::~FastqRecord() {
    if (aln != NULL) {
        bam_destroy1(aln);
        aln = NULL;
    }
}

/* Write a FastqRecord object to a FASTQ file */
int FastqRecord::write(ogzstream& out_fq, std::string name) {
    if (aln != NULL) {
        seq_str = get_forward_read(aln);
        std::string qual_seq("");
        uint8_t* qual = bam_get_qual(aln);
        if (qual[0] == 255)
            qual_seq = "*";
        else {
            for (auto i = 0; i < aln->core.l_qseq; ++i) {
                qual_seq += (char)(qual[i] + 33);
            }
        }
        if (aln->core.flag & BAM_FREVERSE)
            std::reverse(qual_seq.begin(), qual_seq.end());
        qual_str = qual_seq;
    }
    out_fq << "@" << name << "\n";
    out_fq << seq_str << "\n+\n";
    out_fq << qual_str << "\n";
    return 1;
}

/* Split a string on a delimiter
 * From: https://stackoverflow.com/a/64886763
 */
std::vector<std::string> str_to_vector(const std::string str,
                                       const std::string regex_str) {
    std::regex regexz(regex_str);
    std::vector<std::string> vec(
        std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
        std::sregex_token_iterator());
    return vec;
}

std::set<std::string> str_to_set(const std::string str,
                                 const std::string regex_str) {
    std::regex regexz(regex_str);
    std::set<std::string> list(
        std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
        std::sregex_token_iterator());
    return list;
}

/// @brief Returns true if a bam record is reversely lifted.
/// @param aln_orig The original bam record
/// @param aln_lft  The lifted bam record
/// @return true if reverse
bool is_reversedly_lifted(bam1_t* aln_orig, bam1_t* aln_lft) {
    return ((aln_orig->core.flag & BAM_FREVERSE) ^
            (aln_lft->core.flag & BAM_FREVERSE));
}

/// @brief Updates BAM tags specific to Ultima Genomics.
/// @param aln aln The bam record to read and update
/// @param rev Whether to reverse the record
void update_ultima_genomics_tags(bam1_t* aln, bool rev) {
    // TP tag
    uint8_t* tp_ptr = bam_aux_get(aln, "tp");
    if (tp_ptr != NULL) {
        std::vector<int8_t> tp_array;
        for (uint32_t i = 0; i < bam_auxB_len(tp_ptr); i++) {
            tp_array.push_back(bam_auxB2i(tp_ptr, i));
        }
        std::reverse(tp_array.begin(), tp_array.end());
        if (bam_aux_update_array(aln, "tp", 'c', tp_array.size(),
                                 tp_array.data()) != 0) {
            std::cerr << "[W::update_ultima_genomics_tags] Failed to update "
                         "the tp tag for read "
                      << bam_get_qname(aln) << "\n";
        }
    } else {
        std::cerr
            << "[W::update_ultima_genomics_tags] No tp tag is found for read "
            << bam_get_qname(aln) << "\n";
    }
    // T0 tag
    uint8_t* t0_ptr = bam_aux_get(aln, "t0");
    if (t0_ptr != NULL) {
        std::string t0_str = bam_aux2Z(t0_ptr);
        std::reverse(t0_str.begin(), t0_str.end());
        if (bam_aux_update_str(aln, "t0", t0_str.length(), t0_str.c_str()) !=
            0) {
            std::cerr << "[W::update_ultima_genomics_tags] Failed to update "
                         "the t0 tag\n";
        }
    } else {
        std::cerr
            << "[W::update_ultima_genomics_tags] No t0 tag is found for read "
            << bam_get_qname(aln) << "\n";
    }
}

/* Update SEQ and QUAL if in a reversed chain */
int reverse_seq_and_qual(bam1_t* aln) {
    bam1_core_t* c = &(aln->core);
    uint8_t* seq = bam_get_seq(aln);
    uint8_t* qual = bam_get_qual(aln);
    std::vector<int> read_nt16;
    std::vector<uint8_t> qual_tmp;
    for (int n = 0; n < c->l_qseq; n++) {
        read_nt16.push_back(seq_comp_table[bam_seqi(seq, n)]);
        qual_tmp.push_back(qual[n]);
    }
    // Make sure the length of the reversed seq is concordant with the original
    // seq
    if (read_nt16.size() != c->l_qseq) {
        return -1;
    }
    for (int n = 0; n < c->l_qseq; n++) {
        bam_set_seqi(seq, n, read_nt16[c->l_qseq - n - 1]);
        qual[n] = qual_tmp[c->l_qseq - n - 1];
    }

    return 0;
}

std::vector<std::pair<std::string, int32_t>> fai_to_map(std::string fai_fn) {
    std::ifstream fai_fp(fai_fn);
    std::string line, name;
    int32_t length;
    std::vector<std::pair<std::string, int32_t>> lengths;
    while (getline(fai_fp, line)) {
        auto split_line = str_to_vector(line, "\t");
        name = split_line[0];
        length = std::stoi(split_line[1]);
        lengths.push_back(std::make_pair(name, length));
    }
    fai_fp.close();
    return lengths;
}

sam_hdr_t* lengthmap_to_hdr(std::vector<std::pair<std::string, int32_t>> lm,
                            const sam_hdr_t* const hdr_orig) {
    sam_hdr_t* hdr = sam_hdr_dup(hdr_orig);
    // Clear all contig length info in the original header
    sam_hdr_remove_lines(hdr, "SQ", "SN", NULL);
    for (auto& it : lm) {
        std::string name = it.first;
        std::string length = std::to_string(it.second);
        if (sam_hdr_add_line(hdr, "SQ", "SN", name.c_str(), "LN",
                             length.c_str(), NULL) < 0) {
            std::cerr << "[W::utils::lengthmap_to_hdr] error during "
                         "updating BAM header\n";
            std::cerr << name << " " << length << "\n";
        }
    }
    return hdr;
}

sam_hdr_t* fai_to_hdr(std::string fai_fn, const sam_hdr_t* const hdr_orig) {
    sam_hdr_t* hdr = sam_hdr_dup(hdr_orig);
    // Clear all contig length info in the original header
    sam_hdr_remove_lines(hdr, "SQ", "SN", NULL);
    std::ifstream fai_fp(fai_fn);
    std::string line;
    std::string name, length;
    while (getline(fai_fp, line)) {
        auto split_line = str_to_vector(line, "\t");
        name = split_line[0];
        length = split_line[1];
        if (sam_hdr_add_line(hdr, "SQ", "SN", name.c_str(), "LN",
                             length.c_str(), NULL) < 0) {
            std::cerr << "[W::utils::fai_to_hdr] error during "
                         "updating BAM header\n";
            std::cerr << name << " " << length << "\n";
        }
    }
    fai_fp.close();
    return hdr;
}

/**
 * Reduces the size of a BAM object, usually after another data change.
 * Adapted from
 * https://github.com/samtools/htslib/blob/4ff46a6f609fbf886457bbab0f3253446b46a541/sam.c#L429
 *
 * @param b A BAM object.
 * @param desired Number of bytes to trim.
 */
void _realloc_bam_data(bam1_t* b, size_t desired) {
    uint32_t new_m_data;
    uint8_t* new_data;
    new_m_data = desired;
    kroundup32(new_m_data);
    if (new_m_data < desired) {
        errno = ENOMEM;  // Not strictly true but we can't store the size
        throw std::runtime_error(
            "Failed to realloc BAM data - cannot allocate memory");
    }
    if ((bam_get_mempolicy(b) & BAM_USER_OWNS_DATA) == 0) {
        new_data = static_cast<uint8_t*>(realloc(b->data, new_m_data));
    } else {
        if ((new_data = static_cast<uint8_t*>(malloc(new_m_data))) != NULL) {
            if (b->l_data > 0)
                memcpy(new_data, b->data,
                       b->l_data < b->m_data ? b->l_data : b->m_data);
            bam_set_mempolicy(b, bam_get_mempolicy(b) & (~BAM_USER_OWNS_DATA));
        }
    }
    if (!new_data) throw std::runtime_error("Failed to realloc BAM data");
    b->data = new_data;
    b->m_data = new_m_data;
}

// Serialize a `vector<pair<string, int32_t>>` object
size_t serialize_lengthmap(
    std::ofstream& out,
    std::vector<std::pair<std::string, int32_t>> length_map) {
    size_t size = 0;
    size_t nelems = length_map.size();
    out.write(reinterpret_cast<char*>(&nelems), sizeof(nelems));
    size += sizeof(nelems);
    for (auto& x : length_map) {
        size_t str_size = x.first.size();
        out.write(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        out.write(reinterpret_cast<const char*>(x.first.data()), str_size);
        out.write(reinterpret_cast<char*>(&x.second), sizeof(x.second));
        size += sizeof(str_size) + str_size + sizeof(x.second);
    }
    return size;
}

// Load a serialized `vector<pair<string, int32_t>>` object
std::vector<std::pair<std::string, int32_t>> load_lengthmap(std::ifstream& in) {
    size_t map_size;
    in.read(reinterpret_cast<char*>(&map_size), sizeof(map_size));
    std::vector<std::pair<std::string, int32_t>> length_map;
    for (auto i = 0; i < map_size; ++i) {
        size_t str_size;
        in.read(reinterpret_cast<char*>(&str_size), sizeof(str_size));
        std::vector<char> buf(str_size);
        in.read(reinterpret_cast<char*>(buf.data()), str_size);
        std::string key(buf.begin(), buf.end());
        int32_t value;
        in.read(reinterpret_cast<char*>(&value), sizeof(value));
        length_map.push_back(std::pair<std::string, int32_t>(key, value));
    }
    return length_map;
}

std::string make_cmd(int argc, char** argv) {
    std::string cmd("");
    for (auto i = 0; i < argc; ++i) {
        cmd += std::string(argv[i]) + " ";
    }
    return cmd;
}

};  // namespace LevioSamUtils
