/*
 * cigar.cpp
 *
 * Cigar-related utilities.
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "cigar.hpp"

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

namespace Cigar {
/**
 * Updates a cigar vector.
 *
 * It is similar to `cigar.push_back(bam_cigar_gen(len, op))`, and additionally
 * performs CIGAR operator reduction, e.g. (1M2M -> 3M), (1D1I -> 1M)
 *
 * @param cigar The CIGAR object to be updated.
 * @param len The length of the new CIGAR element.
 * @param op The operator of the new CIGAR element.
 * @param no_reduce Set to true to disable CIGAR operator reduction.
 */
void push_cigar(CigarVector& cigar_vec, uint32_t len, uint16_t op,
                const bool no_reduce = false) {
    if (len == 0) return;
    if (cigar_vec.size() == 0 || no_reduce) {
        cigar_vec.push_back(bam_cigar_gen(len, op));
        return;
    }
    auto back_op = bam_cigar_op(cigar_vec.back());
    auto back_type = bam_cigar_type(back_op);
    auto back_len = bam_cigar_oplen(cigar_vec.back());
    auto op_type = bam_cigar_type(op);
    // If operators are the same, merge
    // We also merge S-I/I-S cases
    // We don't merge N-D, H-P because they might mean differently
    if (back_op == op || (back_type == 1 && op_type == 1)) {
        len += back_len;
        cigar_vec.back() = bam_cigar_gen(len, back_op);
        // Cancel out complementary operators
    } else if ((back_type == 2 && op_type == 1) ||
               (back_type == 1 && op_type == 2)) {
        // For cases S-D, we ignore Ds, since deletions don't
        // contribute to either the alignment or the position
        // in soft-clipped regions.
        if (back_op == BAM_CSOFT_CLIP) {
            return;
        }
        if (len == back_len) {
            cigar_vec.pop_back();
            push_cigar(cigar_vec, len, BAM_CMATCH, false);
        } else if (len > back_len) {
            len -= back_len;
            cigar_vec.back() = bam_cigar_gen(back_len, BAM_CMATCH);
            cigar_vec.push_back(bam_cigar_gen(len, op));
        } else {
            back_len -= len;
            cigar_vec.back() = bam_cigar_gen(back_len, back_op);
            cigar_vec.push_back(bam_cigar_gen(len, BAM_CMATCH));
        }
    } else
        cigar_vec.push_back(bam_cigar_gen(len, op));
}

/**
 * Pops `size` bases (from the back) of a cigar vector.
 *
 * @param cigar_vec
 * @param size
 */
void pop_cigar(CigarVector& cigar_vec, uint32_t size) {
    while (size > 0) {
        auto cg = cigar_vec.back();
        cigar_vec.pop_back();
        auto cigar_op_len = bam_cigar_oplen(cg);
        auto cigar_op = bam_cigar_op(cg);
        // Check if the last operator consumes REF; if it doesn't, just pop it
        if (bam_cigar_type(cigar_op) & 1) {
            if (cigar_op_len > size) {
                cigar_op_len -= size;
                size = 0;
                push_cigar(cigar_vec, cigar_op_len, cigar_op, false);
            } else if (cigar_op_len == size) {
                size = 0;
            } else {
                size -= cigar_op_len;
            }
        }
    }
    // TODO: handles the scenario where `size` is greater than
    // the length of the cigar vector.
}

void sclip_cigar_front(uint32_t* cigar, const uint32_t& n_cigar, int len_clip,
                       CigarVector& new_cigar_vec, int& idx,
                       int& query_offset) {
    if (len_clip == 0) return;

    for (auto i = idx; i < n_cigar; i++) {
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        // Only replace bases that consume the QUERY with the SOFT_CLIP ("S")
        // operator
        if (bam_cigar_type(cigar_op) & 1) {
            // If a CIGAR OP is longer than #clipped, append clipped bases to
            // the updated CIGAR and truncate the original CIGAR for future
            // usage.
            if (cigar_op_len >= len_clip) {
                cigar_op_len -= len_clip;
                push_cigar(new_cigar_vec, len_clip, BAM_CSOFT_CLIP, false);
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
                push_cigar(new_cigar_vec, cigar_op_len, BAM_CSOFT_CLIP, false);
                query_offset += cigar_op_len;
                len_clip -= cigar_op_len;
            }
        }
    }
}

void sclip_cigar_back(CigarVector& cigar_vec, int len_clip) {
    if (len_clip == 0) return;

    int remaining_lc = len_clip;
    for (int i = cigar_vec.size() - 1; i >= 0; i--) {
        auto cigar_op_len = bam_cigar_oplen(cigar_vec[i]);
        auto cigar_op = bam_cigar_op(cigar_vec[i]);
        cigar_vec.pop_back();
        // Only replace bases that consume the QUERY with the SOFT_CLIP ("S")
        // operator
        if (bam_cigar_type(cigar_op) & 1) {
            if (cigar_op_len >= remaining_lc) {
                // If a CIGAR OP is longer than #clipped, append clipped bases
                // to the updated CIGAR and truncate the original CIGAR.
                cigar_op_len -= remaining_lc;
                push_cigar(cigar_vec, cigar_op_len, cigar_op, false);
                push_cigar(cigar_vec, len_clip, BAM_CSOFT_CLIP, false);
                return;
            } else if (cigar_op == BAM_CMATCH) {
                // If a CIGAR OP is not long enough, shorten `remaining_lc`
                // and proceed to an earlier OP
                remaining_lc -= cigar_op_len;
            } else if (cigar_op == BAM_CINS || cigar_op == BAM_CSOFT_CLIP) {
                len_clip += cigar_op_len;
            } else {
                std::cerr
                    << "[W::chain::sclip_cigar_back] Unexpected OP during "
                       "back_clipping: "
                    << cigar_op << "\n";
            }
        }
    }
}

/**
 * Sets the cigar string of a BAM record to be empty ("*").
 *
 * @param aln A BAM object.
 */
void set_empty_cigar(bam1_t* aln) {
    uint32_t prev_n_cigar = aln->core.n_cigar;
    size_t new_m_data = (size_t)aln->l_data - prev_n_cigar * 4;

    // Move data to the correct place.
    // Info in the `data` field:
    //   `c = aln->core; d = aln->data;`
    //   qname: `d` to `d + c.l_qname`
    //   cigar: `d + c.l_qname` to `d + c.l_qname + c.n_cigar<<2`
    // Here everything after the old CIGAR is moved to right after QNAME.
    memmove(aln->data + aln->core.l_qname,
            aln->data + aln->core.l_qname + prev_n_cigar * 4,
            new_m_data - aln->core.l_qname);
    _realloc_bam_data(aln, new_m_data);
    aln->core.n_cigar = 0;
    aln->l_data -= prev_n_cigar * 4;
}

void update_cigar(bam1_t* aln, CigarVector& new_cigar_vec) {
    uint32_t* cigar = bam_get_cigar(aln);
    // Adapted from samtools/sam.c
    // https://github.com/samtools/htslib/blob/2264113e5df1946210828e45d29c605915bd3733/sam.c#L515
    if (aln->core.n_cigar != new_cigar_vec.size()) {
        auto cigar_st = (uint8_t*)bam_get_cigar(aln) - aln->data;
        auto fake_bytes = aln->core.n_cigar * 4;
        aln->core.n_cigar = (uint32_t)new_cigar_vec.size();
        auto n_cigar4 = aln->core.n_cigar * 4;
        auto orig_len = aln->l_data;
        if (n_cigar4 > fake_bytes) {
            // Check if we need to update `aln->m_data`.
            //
            // Adapted from htslib/sam_internal.h
            // https://github.com/samtools/htslib/blob/31f0a76d338c9bf3a6893b71dd437ef5fcaaea0e/sam_internal.h#L48
            auto new_m_data = (size_t)aln->l_data + n_cigar4 - fake_bytes;
            kroundup32(new_m_data);
            if (new_m_data > aln->m_data) {
                auto new_data =
                    static_cast<uint8_t*>(realloc(aln->data, new_m_data));
                if (!new_data) {
                    std::cerr << "[Error] Failed to expand a bam1_t struct for "
                              << bam_get_qname(aln) << "\n";
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
        if (n_cigar4 > fake_bytes) {
            memcpy(aln->data + cigar_st,
                   aln->data + (n_cigar4 - fake_bytes) + 8, n_cigar4);
        }
        aln->core.n_cigar = new_cigar_vec.size();
    }
    for (int i = 0; i < aln->core.n_cigar; i++) {
        *(cigar + i) = new_cigar_vec[i];
    }
}

void debug_print_cigar(uint32_t* cigar, size_t n_cigar) {
    for (int i = 0; i < n_cigar; i++) {
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        std::cerr << cigar_op_len << bam_cigar_opchr(cigar_op);
    }
    std::cerr << "\n";
}

}  // namespace Cigar