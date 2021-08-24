#ifndef LEVIOSAM_UTILS_HPP
#define LEVIOSAM_UTILS_HPP

#include <iostream>
#include <fstream>
#include <htslib/sam.h>
#include <vector>

const int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 };

namespace LevioSamUtils {
class WriteToFastq {
public:
    WriteToFastq() {};
    ~WriteToFastq() {
        out_fqS.close();
        out_fq1.close();
        out_fq2.close();
    };

    void init(
        const std::string outpre, const std::string sm,
        const int mc
    );

    std::ofstream out_fqS;
    std::ofstream out_fq1;
    std::ofstream out_fq2;
    std::string split_mode = "";
    std::mutex mutex_fwrite_fq;
    int mapq_cutoff;
};

void update_cigar(bam1_t* aln, std::vector<uint32_t> &new_cigar);
void debug_print_cigar(uint32_t* cigar, size_t n_cigar);
void remove_mn_md_tag(bam1_t* aln);
static std::string get_read(const bam1_t *rec);
void write_fq_from_bam_core(bam1_t* aln, std::ofstream& out_fq);
void write_fq_from_bam(bam1_t* aln, WriteToFastq* w2fq);
}

#endif
