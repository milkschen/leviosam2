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
    const std::string& deferred_sam_fname,
    const std::string& out_deferred_sam_fname);

class WriteToFastq {
public:
    WriteToFastq() {
        write_deferred = false;
    };
    ~WriteToFastq();

    void init(
        const std::string outpre, const std::string sm,
        const int mapq, const int isize,
        const float clipped_frac, const int aln_score,
        const std::string of
    );

    void write_fq_from_bam(bam1_t* aln);
    void write_low_mapq_bam(bam1_t* aln, bam_hdr_t* hdr);

    std::ofstream out_fqS;
    std::ofstream out_fq1;
    std::ofstream out_fq2;
    samFile* out_fp;
    std::string split_mode = "";
    std::mutex mutex_fwrite_fq;
    int min_mapq;
    int max_isize;
    float max_clipped_frac;
    int min_aln_score;
    std::unordered_map<std::string, FastqRecord> r1_db;
    std::unordered_map<std::string, FastqRecord> r2_db;
    bool write_deferred;
};

void update_cigar(bam1_t* aln, std::vector<uint32_t> &new_cigar);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);
void remove_mn_md_tag(bam1_t* aln);
static std::string get_read(const bam1_t *rec);
void write_fq_from_bam_core(bam1_t* aln, std::ofstream& out_fq);

std::vector<std::string> split_str(
    const std::string str, const std::string regex_str);
}

#endif
