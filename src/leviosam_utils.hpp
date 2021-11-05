#ifndef LEVIOSAM_UTILS_HPP
#define LEVIOSAM_UTILS_HPP

#include <iostream>
#include <fstream>
#include <mutex>
#include <htslib/sam.h>
#include <unordered_map>
#include <vector>
#include <cstring>
#include "robin_hood.h"

const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

namespace LevioSamUtils {

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

    int write(std::ofstream& out_fq, std::string name);
    std::string seq_str;
    std::string qual_str;
    bam1_t* aln = NULL;
};

robin_hood::unordered_map<std::string, FastqRecord> read_unpaired_fq(
    const std::string& fq_fname);
robin_hood::unordered_map<std::string, FastqRecord> read_deferred_bam(
    samFile* dsam_fp, samFile* out_dsam_fp, bam_hdr_t* hdr,
    std::ofstream& out_r1_fp, std::ofstream& out_r2_fp);

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
        const std::string of, bam_hdr_t* ihdr
    );

    void write_deferred_bam(bam1_t* aln, bam_hdr_t* hdr);
    void write_deferred_bam_orig(bam1_t* aln);

    samFile* out_fp;
    samFile* out_fp_orig;
    bam_hdr_t* hdr_orig;
    std::string split_mode = "";
    std::mutex mutex_fwrite;
    int min_mapq;
    int max_isize;
    float max_clipped_frac;
    int min_aln_score;
    bool write_deferred;
};

void update_cigar(bam1_t* aln, std::vector<uint32_t> &new_cigar);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);
void remove_mn_md_tag(bam1_t* aln);
static std::string get_read(const bam1_t *rec);

std::vector<std::string> split_str(
    const std::string str, const std::string regex_str);
int reverse_seq_and_qual(bam1_t* aln);
sam_hdr_t* fai_to_hdr(std::string dest_fai_fname, const sam_hdr_t* const hdr_orig);

}

#endif
