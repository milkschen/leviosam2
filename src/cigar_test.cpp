/*
 * ciar_test.cpp
 *
 * Test cigar related utilities for levioSAM2
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "cigar.hpp"

#include <iostream>

#include "gtest/gtest.h"
#include "leviosam_utils.hpp"

TEST(CigarTest, PushCigar) {
    Cigar::CigarVector cigar1;
    Cigar::push_cigar(cigar1, 2, BAM_CMATCH, false);
    Cigar::push_cigar(cigar1, 1, BAM_CMATCH, false);
    EXPECT_EQ(cigar1.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar1.back()), 3);
    EXPECT_EQ(bam_cigar_op(cigar1.back()), BAM_CMATCH);

    Cigar::CigarVector cigar2;
    Cigar::push_cigar(cigar2, 1, BAM_CINS, false);
    Cigar::push_cigar(cigar2, 1, BAM_CDEL, false);
    EXPECT_EQ(cigar2.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar2.back()), 1);
    EXPECT_EQ(bam_cigar_op(cigar2.back()), BAM_CMATCH);
}

TEST(CigarTest, SClipCigarFront1) {
    // 5M -> 3S2M
    Cigar::CigarVector cigar;
    Cigar::push_cigar(cigar, 5, BAM_CMATCH, false);
    Cigar::CigarVector new_cigar;
    int idx = 0, q = 0;
    Cigar::sclip_cigar_front(&cigar[0], cigar.size(), 3, new_cigar, idx, q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 3);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 0);
    EXPECT_EQ(q, 3);
}

TEST(CigarTest, SClipCigarFront2) {
    // 4S2M -> 5S1M
    Cigar::CigarVector cigar;
    Cigar::push_cigar(cigar, 4, BAM_CSOFT_CLIP, false);
    Cigar::push_cigar(cigar, 2, BAM_CMATCH, false);
    Cigar::CigarVector new_cigar;
    int idx = 0, q = 0;
    Cigar::sclip_cigar_front(&cigar[0], cigar.size(), 5, new_cigar, idx, q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 5);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 1);
    EXPECT_EQ(q, 5);
}

TEST(CigarTest, SClipCigarBack) {
    Cigar::CigarVector cigar1;
    Cigar::push_cigar(cigar1, 5, BAM_CMATCH, false);
    Cigar::CigarVector new_cigar1;
}

TEST(CigarTest, SetEmptyCigar) {
    std::string hdr_str = "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:10000";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();

    std::string record1 =
        "read1\t81\tchr1\t100\t33\t6S10M\t"
        "=\t300\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:10M\tNM:i:0";
    kstring_t str1;
    str1.s = (char*)record1.c_str();
    str1.l = record1.length();
    str1.m = kstr_get_m(str1.l);
    EXPECT_GE(sam_parse1(&str1, sam_hdr, aln), 0);

    Cigar::set_empty_cigar(aln);
    std::string seq = "ATTACATTCCATTCCA";
    EXPECT_EQ(aln->core.n_cigar, 0);
    // SEQ is after CIGAR, so lets check if SEQ looks good.
    for (auto i = 0; i < aln->core.l_qseq; i++) {
        EXPECT_EQ(seq_nt16_str[bam_seqi(bam_get_seq(aln), i)], seq[i]);
    }
    bam_destroy1(aln);

    // Second test case. The original CIGAR is more complicated.
    aln = bam_init1();
    record1 =
        "read2\t81\tchr1\t100\t33\t2S3M2I2M2D7M1H\t"
        "=\t300\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:10M\tNM:i:0";
    str1.s = (char*)record1.c_str();
    str1.l = record1.length();
    str1.m = kstr_get_m(str1.l);
    EXPECT_EQ(sam_parse1(&str1, sam_hdr, aln), 0);
    Cigar::set_empty_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 0);
    for (auto i = 0; i < aln->core.l_qseq; i++) {
        EXPECT_EQ(seq_nt16_str[bam_seqi(bam_get_seq(aln), i)], seq[i]);
    }
    bam_destroy1(aln);
}

TEST(CigarTest, UpdateCigar) {
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();

    kstring_t str1;
    std::string record1 =
        "read1\t81\tchr1\t145334831\t33\t6S10M\t"
        "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:10M";
    str1.s = (char*)record1.c_str();
    str1.l = record1.length();
    str1.m = kstr_get_m(str1.l);
    EXPECT_EQ(sam_parse1(&str1, sam_hdr, aln), 0);

    Cigar::CigarVector new_cigar{bam_cigar_gen(16, BAM_CMATCH)};
    Cigar::update_cigar(aln, new_cigar);
    EXPECT_EQ(aln->core.n_cigar, 1);
    EXPECT_EQ(*bam_get_cigar(aln), bam_cigar_gen(16, BAM_CMATCH));
    bam_destroy1(aln);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
