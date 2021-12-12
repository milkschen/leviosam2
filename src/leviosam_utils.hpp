/*
 * leviosam_utils.hpp
 *
 * Utility functions for the leviosam software
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#ifndef LEVIOSAM_UTILS_HPP
#define LEVIOSAM_UTILS_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <htslib/sam.h>
#include <set>
#include <vector>
#include <cstring>
#include "bed.hpp"
#include "gzstream.h"
#include "robin_hood.h"

const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

namespace LevioSamUtils {

const std::vector<std::string> DEFER_OPT
    {"lifted", "mapq", "clipped_frac", "isize", "aln_score", "hdist"};

class FastqRecord {
public:
    ~FastqRecord();
    FastqRecord(bam1_t* a);
    FastqRecord(const std::string& seq, const std::string& qual);

    // Copy constructor
    FastqRecord(const FastqRecord& rhs):
        seq_str(rhs.seq_str), qual_str(rhs.qual_str) {
        if (rhs.aln != NULL) {
            aln = bam_init1();
            if (bam_copy1(aln, rhs.aln) == NULL) {
                std::cerr << "Warning: copying failed for " << bam_get_qname(rhs.aln) << "\n";
            }
        }
    }

    // Move constructor
    FastqRecord(FastqRecord&& rhs) noexcept:
        seq_str(rhs.seq_str), qual_str(rhs.qual_str) {
        if (rhs.aln != NULL) {
            aln = bam_init1();
            if (bam_copy1(aln, rhs.aln) == NULL) {
                std::cerr << "Warning: moving failed for " << bam_get_qname(rhs.aln) << "\n";
                bam_destroy1(rhs.aln);
                rhs.aln = NULL;
            }
        }
    }

    int write(ogzstream& out_fq, std::string name);
    std::string seq_str;
    std::string qual_str;
    bam1_t* aln = NULL;
};

typedef robin_hood::unordered_map<std::string, FastqRecord> fastq_map;

fastq_map read_unpaired_fq(
    const std::string& fq_fname);
fastq_map read_deferred_bam(
    samFile* dsam_fp, samFile* out_dsam_fp, sam_hdr_t* hdr,
    ogzstream& out_r1_fp, ogzstream& out_r2_fp);

class WriteDeferred {
public:
    WriteDeferred() {
        write_deferred = false;
    };
    ~WriteDeferred();

    void init(
        const std::string outpre, const std::string sm,
        const int mapq, const int isize,
        const float clipped_frac, const int aln_score,
        const std::string of, sam_hdr_t* ihdr, sam_hdr_t* ohdr,
        const BedUtils::Bed &b_defer_source,
        const BedUtils::Bed &b_defer_dest,
        const BedUtils::Bed &b_commit_source,
        const BedUtils::Bed &b_commit_dest
    );

    void print_info();
    void write_deferred_bam(bam1_t* aln, sam_hdr_t* hdr);
    void write_deferred_bam_orig(bam1_t* aln);
    bool commit_aln_dest(const bam1_t* const aln);
    bool commit_aln_source(const bam1_t* const aln);

    std::mutex mutex_fwrite;

private:
    bool write_deferred;
    samFile* out_fp;
    samFile* out_fp_orig;
    sam_hdr_t* hdr;
    sam_hdr_t* hdr_orig;
    std::set<std::string> split_modes;
    int min_mapq, max_isize, min_aln_score, max_hdist;
    float max_clipped_frac;
    BedUtils::Bed bed_defer_source, bed_defer_dest, bed_commit_source, bed_commit_dest;
};

void update_cigar(bam1_t* aln, std::vector<uint32_t> &new_cigar);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);
void remove_mn_md_tag(bam1_t* aln);
static std::string get_read(const bam1_t *rec);

std::vector<std::string> str_to_vector(
    const std::string str, const std::string regex_str);
std::set<std::string> str_to_set(
    const std::string str, const std::string regex_str);

int reverse_seq_and_qual(bam1_t* aln);
sam_hdr_t* fai_to_hdr(std::string fai_fn, const sam_hdr_t* const hdr_orig);
std::vector<std::pair<std::string, int32_t>> fai_to_map(std::string fai_fn);
sam_hdr_t* lengthmap_to_hdr(
    std::vector<std::pair<std::string, int32_t>> lm,
    const sam_hdr_t* const hdr_orig);

std::vector<std::pair<std::string, int32_t>> load_lengthmap(std::ifstream& in);
size_t serialize_lengthmap(
    std::ofstream& out,
    std::vector<std::pair<std::string, int32_t>> length_map);

std::string make_cmd(int argc, char** argv);

}

#endif
