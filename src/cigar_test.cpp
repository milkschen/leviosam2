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

// #include <unistd.h>

#include <iostream>

#include "gtest/gtest.h"
// #include "leviosam_utils.hpp"

TEST(CigarTest, PushCigar) {
    std::vector<uint32_t> cigar1;
    Cigar::push_cigar(cigar1, 2, BAM_CMATCH, false);
    Cigar::push_cigar(cigar1, 1, BAM_CMATCH, false);
    EXPECT_EQ(cigar1.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar1.back()), 3);
    EXPECT_EQ(bam_cigar_op(cigar1.back()), BAM_CMATCH);

    std::vector<uint32_t> cigar2;
    Cigar::push_cigar(cigar2, 1, BAM_CINS, false);
    Cigar::push_cigar(cigar2, 1, BAM_CDEL, false);
    EXPECT_EQ(cigar2.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar2.back()), 1);
    EXPECT_EQ(bam_cigar_op(cigar2.back()), BAM_CMATCH);
}

TEST(CigarTest, SClipCigarFront1) {
    // 5M -> 3S2M
    std::vector<uint32_t> cigar;
    Cigar::push_cigar(cigar, 5, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar;
    int idx = 0, q = 0;
    Cigar::sclip_cigar_front(&cigar[0], cigar.size(), 3, new_cigar, idx,
                                     q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 3);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 0);
    EXPECT_EQ(q, 3);
}

TEST(CigarTest, SClipCigarFront2) {
    // 4S2M -> 5S1M
    std::vector<uint32_t> cigar;
    Cigar::push_cigar(cigar, 4, BAM_CSOFT_CLIP, false);
    Cigar::push_cigar(cigar, 2, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar;
    int idx = 0, q = 0;
    Cigar::sclip_cigar_front(&cigar[0], cigar.size(), 5, new_cigar, idx,
                                     q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 5);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 1);
    EXPECT_EQ(q, 5);
}

TEST(CigarTest, SClipCigarBack) {
    std::vector<uint32_t> cigar1;
    Cigar::push_cigar(cigar1, 5, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar1;
}

// TEST(UtilsTest, UpdateCigar) {
//     std::string hdr_str =
//         "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
//     sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
//     bam1_t* aln = bam_init1();
//     int err;

//     kstring_t str1;
//     std::string record1 =
//         "read1\t81\tchr1\t145334831\t33\t6S10M\t"
//         "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:10M";
//     str1.s = (char*)record1.c_str();
//     str1.l = record1.length();
//     str1.m = kstr_get_m(str1.l);
//     err = sam_parse1(&str1, sam_hdr, aln);
//     EXPECT_EQ(err, 0);

//     std::vector<uint32_t> new_cigar;
//     for (int i = 0; i < aln->core.n_cigar; i++) {
//         auto cigar_op_len = bam_cigar_oplen(cigar[i]);
//         auto cigar_op = bam_cigar_op(cigar[i]);
//         push_cigar(new_cigar, cigar_op_len, cigar_op, false);
//     }
//     EXPECT_EQ(LevioSamUtils::update_cigar(aln), "16M");
//     bam_destroy1(aln);
// }

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
