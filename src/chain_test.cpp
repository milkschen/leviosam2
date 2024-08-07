/*
 * chain_test.cpp
 *
 * Test chain related functions for levioSAM2
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "chain.hpp"

#include <unistd.h>

#include <iostream>

#include "gtest/gtest.h"
#include "leviosam_utils.hpp"

/* Chain tests */
TEST(ChainTest, SimpleRankAndLift) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);

    int pos = 674047;
    std::string contig = "chr1";
    int rank = cmap.get_start_rank(contig, pos);
    EXPECT_EQ(rank, 1);
    EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
    EXPECT_EQ(cmap.lift_pos(contig, pos, 0, true), 100272);

    pos = 207130;
    contig = "chr2";
    rank = cmap.get_start_rank(contig, pos);
    EXPECT_EQ(rank, 2);
    EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
    EXPECT_EQ(cmap.lift_pos(contig, pos, 0, true), 206846 + 121 + 8);
}

TEST(ChainTest, ParseChainLineCornerZero) {
    std::string hdr =
        "chain 5 corner_zero 28 + 0 28 corner_zero_dest 300 + 100 130 0";
    std::string line1 = "20\t0\t2\n";
    std::string line2 = "8\n";
    chain::BitVectorMap start_bv_map, end_bv_map;
    std::string source, target;
    int source_offset = 0, target_offset = 0, source_len = 0;
    bool current_ss = true;
    chain::ChainMap cmap;
    chain::LengthMap lmap{std::make_pair("corner_zero_dest", 300)};
    cmap.parse_chain_line(hdr, source, target, source_len, source_offset,
                          target_offset, current_ss, start_bv_map, end_bv_map,
                          lmap);
    cmap.parse_chain_line(line1, source, target, source_len, source_offset,
                          target_offset, current_ss, start_bv_map, end_bv_map,
                          lmap);
    cmap.parse_chain_line(line2, source, target, source_len, source_offset,
                          target_offset, current_ss, start_bv_map, end_bv_map,
                          lmap);

    for (auto &i : start_bv_map) {
        EXPECT_EQ(i.first, "corner_zero");
    }
    for (int i = 0; i < 28; i++) {
        if (i == 0 || i == 19)
            EXPECT_EQ(start_bv_map["corner_zero"][i], 1);
        else
            EXPECT_EQ(start_bv_map["corner_zero"][i], 0);
        if (i == 20 || i == 27)
            EXPECT_EQ(end_bv_map["corner_zero"][i], 1);
        else
            EXPECT_EQ(end_bv_map["corner_zero"][i], 0);
    }
}

TEST(ChainTest, LiftInReversedRegion) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387497));
    chain::ChainMap cmap("chr1_reversed_region.chain", 0, 0, lm);

    std::string contig = "chr1";
    int pos_array[4] = {146735453, 146735605, 146735135, 146735235};
    int gold_pos_array[4] = {148073114, 148072962, 148073432, 148073332};
    int gold_rank = 154;
    int rank;
    int pos;
    for (int i = 0; i < 4; i++) {
        pos = pos_array[i];
        rank = cmap.get_start_rank(contig, pos);
        EXPECT_EQ(rank, gold_rank);
        EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
        EXPECT_EQ(cmap.lift_pos(contig, pos, 0, true), gold_pos_array[i]);
    }
}

TEST(ChainTest, LiftBamInReversedRegion) {
    samFile *gold_fp = sam_open("HG002-0.3x-bwa-grch37-chr1_rev.sam", "r");
    sam_hdr_t *hdr_gold = sam_hdr_read(gold_fp);
    bam1_t *aln_gold = bam_init1();
    std::vector<bam1_t *> gold1, gold2;
    while (sam_read1(gold_fp, hdr_gold, aln_gold) == 0) {
        if (aln_gold->core.flag & BAM_FREAD1)
            gold1.push_back(bam_dup1(aln_gold));
        else
            gold2.push_back(bam_dup1(aln_gold));
    }
    EXPECT_EQ(gold1.size(), 2);
    EXPECT_EQ(gold2.size(), 2);

    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 249250621));
    chain::ChainMap cmap("hg38_to_hg19-chr1_100216456_104974084.chain", 0, 0,
                         lm);
    samFile *sam_fp = sam_open("HG002-0.3x-bwa-grch38-chr1_rev.sam", "r");
    sam_hdr_t *hdr = sam_hdr_read(sam_fp);
    bam1_t *aln = bam_init1();
    int i = 0;
    while (sam_read1(sam_fp, hdr, aln) == 0) {
        std::string qname = bam_get_qname(aln);
        std::string dest_contig = hdr->target_name[aln->core.tid];
        EXPECT_EQ(cmap.lift_segment(aln, hdr, hdr_gold, true, dest_contig),
                  true);
        auto pos = aln->core.pos;
        if (aln->core.flag & BAM_FREAD1) {
            EXPECT_EQ(qname, bam_get_qname(gold1[i / 2]));
            EXPECT_EQ(pos, gold1[i / 2]->core.pos);
        } else {
            EXPECT_EQ(qname, bam_get_qname(gold2[i / 2]));
            EXPECT_EQ(pos, gold2[i / 2]->core.pos);
        }
        i++;
    }
}

TEST(ChainTest, LiftCigarCoreOneRun) {
    std::vector<uint32_t> new_cigar;
    int query_offset = 0;
    int tmp_gap = 0;
    uint32_t qlen = 252;
    chain::ChainMap cmap;
    std::queue<std::tuple<int32_t, int32_t>> bp;
    bp.push(std::make_tuple(195, -18));
    bp.push(std::make_tuple(217, -2));

    cmap.lift_cigar_core_one_run(new_cigar, bp, 200, BAM_CMATCH, qlen, tmp_gap,
                                 query_offset);
    EXPECT_EQ(tmp_gap, 13);
    cmap.lift_cigar_core_one_run(new_cigar, bp, 5, BAM_CINS, qlen, tmp_gap,
                                 query_offset);
    EXPECT_EQ(tmp_gap, 15);
    cmap.lift_cigar_core_one_run(new_cigar, bp, 47, BAM_CMATCH, qlen, tmp_gap,
                                 query_offset);
    EXPECT_EQ(tmp_gap, 0);
    // LevioSamUtils::debug_print_cigar(new_cigar.data(), new_cigar.size());
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 195);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CMATCH);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[1]), 25);
    EXPECT_EQ(bam_cigar_op(new_cigar[1]), BAM_CINS);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[2]), 32);
    EXPECT_EQ(bam_cigar_op(new_cigar[2]), BAM_CMATCH);
}

TEST(ChainTest, LiftCigar1) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // CIGAR should not change
    kstring_t str;
    std::string record =
        "unchanged\t0\tchr1\t674850\t42\t7M13D13M\t"
        "*\t0\t0\tCAGTTTGTAGTATCTGCAAG\t~~~~~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(13, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(13, BAM_CMATCH));
}

TEST(ChainTest, LiftCigar2) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Add 3 BAM_CSOFT_CLIPs in the beginning
    kstring_t str;
    std::string record =
        "16M_3S13M\t0\tchr1\t687455\t42\t16M\t*\t"
        "0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(13, BAM_CMATCH));
}

TEST(ChainTest, LiftCigar3) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Add 2 BAM_CINSs in the middle
    kstring_t str;
    std::string record =
        "16M_10M2I4M\t0\tchr1\t674141\t42\t16M\t*\t"
        "0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(2, BAM_CINS));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(4, BAM_CMATCH));
}

TEST(ChainTest, LiftCigar4) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Add 3 BAM_CDELs in the middle
    kstring_t str;
    std::string record =
        "16M_10M3D6M\t0\tchr1\t674820\t42\t16M\t*\t"
        "0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(3, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(6, BAM_CMATCH));
}

TEST(ChainTest, LiftCigar5) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Add 3 BAM_CSOFT_CLIPs in the beginning (secondary alignment)
    kstring_t str;
    std::string record =
        "16M_3S13M\t256\tchr1\t687455\t42\t"
        "16M\t*\t0\t0\t*\t*";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(13, BAM_CMATCH));
}

TEST(ChainTest, LiftExtendedCigar1) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Replace extended CIGAR ops (`X` and `=`) with `M`
    // CIGAR should not change
    kstring_t str;
    std::string record =
        "7=13D6=1X6=_7M13D13M\t0\tchr1\t674850\t"
        "42\t7=13D6=1X6=\t*\t0\t0\t"
        "CAGTTTGTAGTATCTGCAAG\t~~~~~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(13, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(13, BAM_CMATCH));
}

TEST(ChainTest, LiftExtendedCigar2) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Add 3 BAM_CSOFT_CLIPs in the beginning
    kstring_t str;
    std::string record =
        "10=1X5=_3S13M\t0\tchr1\t687455\t42\t10=1X5=\t*\t"
        "0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(13, BAM_CMATCH));
}

TEST(ChainTest, LiftExtendedCigarReverse1) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387497));
    chain::ChainMap cmap("chr1_reversed_region.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Replace extended CIGAR ops (`X` and `=`) with `M`
    // CIGAR should not change
    kstring_t str;
    std::string record =
        "10=1X5=_16M\t0\tchr1\t145302531\t42\t10=1X5=\t"
        "*\t0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 1);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(16, BAM_CMATCH));
}

TEST(ChainTest, LiftExtendedCigarReverse2) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387497));
    chain::ChainMap cmap("chr1_reversed_region.chain", 0, 0, lm);
    // chain::ChainMap cmap ("chr1_reversed_region.chain", 5, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Replace extended CIGAR ops (`X` and `=`) with `M`
    // Reverse & 1-bp DEL
    kstring_t str;
    std::string record =
        "10=1X5=_10M1D6M\t0\tchr1\t145331505\t42\t"
        "10=1X5=\t"
        "*\t0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(1, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(6, BAM_CMATCH));
}

TEST(ChainTest, LiftExtendedCigarReverse3) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387497));
    chain::ChainMap cmap("chr1_reversed_region.chain", 0, 0, lm);
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    sam_hdr_t *sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t *aln = bam_init1();
    int err;
    size_t x;
    uint32_t *test_cigar;

    // Replace extended CIGAR ops (`X` and `=`) with `M`
    // Replacing 6-bp MATCH with INS
    kstring_t str;
    std::string record =
        "10=1X5=_7M6I3M\t0\tchr1\t145334831\t42\t10=1X5=\t"
        "*\t0\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~";
    str.s = (char *)record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    // Note: can use the helper function to print out CIGAR results
    // LevioSamUtils::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(6, BAM_CINS));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(3, BAM_CMATCH));
}

TEST(ChainTest, CheckMultiIntvlLegality) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    int sidx = 0;
    int eidx = 1;
    EXPECT_EQ(cmap.check_multi_intvl_legality("chr1", "read1", sidx, eidx, 0),
              false);
    EXPECT_EQ(cmap.check_multi_intvl_legality("chr1", "read1", sidx, eidx, 1),
              false);
    EXPECT_EQ(cmap.check_multi_intvl_legality("chr1", "read1", sidx, eidx, 2),
              true);
    EXPECT_EQ(cmap.check_multi_intvl_legality("chr1", "read1", sidx, eidx, 3),
              true);
}

TEST(ChainTest, UpdateIntervalIndexes) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap("small.chain", 0, 0, lm);
    // cmap.debug_print_intervals("chr1", 5);
    int sidx = 0, eidx = 0;

    // Straightforward test case.
    EXPECT_EQ(cmap.update_interval_indexes("chr1", 674144, sidx, eidx), true);
    EXPECT_EQ(sidx, 0);
    EXPECT_EQ(eidx, 0);

    // Contig not in the map.
    EXPECT_EQ(cmap.update_interval_indexes("chr100", 0, sidx, eidx), false);
    EXPECT_EQ(sidx, -1);
    EXPECT_EQ(eidx, -1);

    // Locus not covered by the chains.
    EXPECT_EQ(cmap.update_interval_indexes("chr1", 674040, sidx, eidx), false);
    EXPECT_EQ(sidx, -1);
    EXPECT_EQ(eidx, 0);

    // Invalid locus: outside chrom
    EXPECT_EQ(cmap.update_interval_indexes("chr1", 1000000000, sidx, eidx),
              false);
    EXPECT_EQ(sidx, -1);
    EXPECT_EQ(eidx, -1);

    EXPECT_EQ(cmap.update_interval_indexes("chr1", 2680090, sidx, eidx), true);
    EXPECT_EQ(sidx, 784);  // TODO: why not 785
    EXPECT_EQ(eidx, 784);  // TODO: why not 785
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
