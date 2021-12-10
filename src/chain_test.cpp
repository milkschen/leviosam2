#include <iostream>
#include <unistd.h>
#include "chain.hpp"
#include "gtest/gtest.h"


size_t kstr_get_m(size_t var) {
    size_t lvar = (size_t)exp2(ceil(log2(var)));
    return lvar;
}

/* Chain tests */
TEST(ChainTest, SimpleRankAndLift) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);

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
    EXPECT_EQ(cmap.lift_pos(contig, pos, 0, true), 206846+121+8);
}


TEST(ChainTest, ParseChainLineCornerZero) {
    std::string hdr = "chain 5 corner_zero 28 + 0 28 corner_zero_dest 300 + 100 130 0";
    std::string line1 = "20\t0\t2\n";
    std::string line2 = "8\n";
    chain::BitVectorMap start_bv_map, end_bv_map;
    std::string source, target;
    int source_offset = 0, target_offset = 0, source_len = 0;
    bool current_ss = true;
    chain::ChainMap cmap;
    cmap.parse_chain_line(
        hdr, source, target, source_len, source_offset, 
        target_offset, current_ss, start_bv_map, end_bv_map);
    cmap.parse_chain_line(
        line1, source, target, source_len, source_offset, 
        target_offset, current_ss, start_bv_map, end_bv_map);
    cmap.parse_chain_line(
        line2, source, target, source_len, source_offset, 
        target_offset, current_ss, start_bv_map, end_bv_map);

    for (auto& i: start_bv_map) {
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
    chain::ChainMap cmap ("chr1_reversed_region.chain", 0, 0, lm);

    std::string contig = "chr1";
    int pos_array [4] = {146735453, 146735605, 146735135, 146735235};
    int gold_pos_array [4] = {148073114, 148072962, 148073432, 148073332};
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
    samFile* gold_fp = sam_open("HG002-0.3x-bwa-grch37-chr1_rev.sam", "r");
    sam_hdr_t* hdr_gold = sam_hdr_read(gold_fp);
    bam1_t* aln_gold = bam_init1();
    std::vector<bam1_t*> gold1, gold2;
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
    chain::ChainMap cmap ("hg38_to_hg19-chr1_100216456_104974084.chain", 0, 0, lm);
    samFile* sam_fp = sam_open("HG002-0.3x-bwa-grch38-chr1_rev.sam", "r");
    sam_hdr_t* hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int i = 0;
    while (sam_read1(sam_fp, hdr, aln) == 0) {
        std::string qname = bam_get_qname(aln);
        std::string dest_contig = hdr->target_name[aln->core.tid];
        EXPECT_EQ(cmap.lift_segment(
            aln, hdr, hdr_gold, true, dest_contig), true);
        auto pos = aln->core.pos;
        if (aln->core.flag & BAM_FREAD1) {
            EXPECT_EQ(qname, bam_get_qname(gold1[i/2]));
            EXPECT_EQ(pos, gold1[i/2]->core.pos);
        } else {
            EXPECT_EQ(qname, bam_get_qname(gold2[i/2]));
            EXPECT_EQ(pos, gold2[i/2]->core.pos);
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

    cmap.lift_cigar_core_one_run(
        new_cigar, bp, 200, BAM_CMATCH,
        qlen, tmp_gap, query_offset);
    EXPECT_EQ(tmp_gap, 13);
    cmap.lift_cigar_core_one_run(
        new_cigar, bp, 5, BAM_CINS,
        qlen, tmp_gap, query_offset);
    EXPECT_EQ(tmp_gap, 15);
    cmap.lift_cigar_core_one_run(
        new_cigar, bp, 47, BAM_CMATCH,
        qlen, tmp_gap, query_offset);
    EXPECT_EQ(tmp_gap, 0);
    // LevioSamUtils::debug_print_cigar(new_cigar.data(), new_cigar.size());
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 195);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CMATCH);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[1]), 25);
    EXPECT_EQ(bam_cigar_op(new_cigar[1]), BAM_CINS);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[2]), 32);
    EXPECT_EQ(bam_cigar_op(new_cigar[2]), BAM_CMATCH);
}


TEST(ChainTest, PushCigar) {
    std::vector<uint32_t> cigar1;
    chain::push_cigar(cigar1, 2, BAM_CMATCH, false);
    chain::push_cigar(cigar1, 1, BAM_CMATCH, false);
    EXPECT_EQ(cigar1.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar1.back()), 3);
    EXPECT_EQ(bam_cigar_op(cigar1.back()), BAM_CMATCH);

    std::vector<uint32_t> cigar2;
    chain::push_cigar(cigar2, 1, BAM_CINS, false);
    chain::push_cigar(cigar2, 1, BAM_CDEL, false);
    EXPECT_EQ(cigar2.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(cigar2.back()), 1);
    EXPECT_EQ(bam_cigar_op(cigar2.back()), BAM_CMATCH);
}


TEST(ChainTest, SClipCigarFront1) {
    // 5M -> 3S2M
    std::vector<uint32_t> cigar;
    chain::push_cigar(cigar, 5, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar;
    int idx = 0, q = 0;
    chain::sclip_cigar_front(&cigar[0], cigar.size(), 3, new_cigar, idx, q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 3);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 0);
    EXPECT_EQ(q, 3);
}


TEST(ChainTest, SClipCigarFront2) {
    // 4S2M -> 5S1M
    std::vector<uint32_t> cigar;
    chain::push_cigar(cigar, 4, BAM_CSOFT_CLIP, false);
    chain::push_cigar(cigar, 2, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar;
    int idx = 0, q = 0;
    chain::sclip_cigar_front(&cigar[0], cigar.size(), 5, new_cigar, idx, q);
    EXPECT_EQ(new_cigar.size(), 1);
    EXPECT_EQ(bam_cigar_oplen(new_cigar[0]), 5);
    EXPECT_EQ(bam_cigar_op(new_cigar[0]), BAM_CSOFT_CLIP);
    EXPECT_EQ(idx, 1);
    EXPECT_EQ(q, 5);
}


TEST(ChainTest, SClipCigarBack) {
    std::vector<uint32_t> cigar1;
    chain::push_cigar(cigar1, 5, BAM_CMATCH, false);
    std::vector<uint32_t> new_cigar1;
}


TEST(ChainTest, LiftCigar1) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    std::string hdr_str = "@HD	VN:1.0	SO:unsorted\n@SQ	SN:chr1	LN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;

    // CIGAR should not change
    kstring_t str;
    std::string record = "unchanged	0	chr1	674850	42	7M13D13M	*	0	0	CAGTTTGTAGTATCTGCAAG	~~~~~~~~~~~~~~~~~~~~";
    str.s = (char*) record.c_str();
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 13, BAM_CMATCH));
}


TEST(ChainTest, LiftCigar2) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    std::string hdr_str = "@HD	VN:1.0	SO:unsorted\n@SQ	SN:chr1	LN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;

    // Add 3 BAM_CSOFT_CLIPs in the beginning
    kstring_t str;
    std::string record = "16M_3S13M	0	chr1	687455	42	16M	*	0	0	ATTACATTCCATTCCA	~~~~~~~~~~~~~~~~";
    str.s = (char*) record.c_str();
    str.l = record.length();
    str.m = kstr_get_m(str.l);
    err = sam_parse1(&str, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    err = cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    EXPECT_EQ(err, 0);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CMATCH));
}


TEST(ChainTest, LiftCigar3) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    std::string hdr_str = "@HD	VN:1.0	SO:unsorted\n@SQ	SN:chr1	LN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;

    // Add 2 BAM_CINSs in the middle
    kstring_t str;
    std::string record = "16M_10M2I4M	0	chr1	674141	42	16M	*	0	0	ATTACATTCCATTCCA	~~~~~~~~~~~~~~~~";
    str.s = (char*) record.c_str();
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 2, BAM_CINS));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 4, BAM_CMATCH));
}


TEST(ChainTest, LiftCigar4) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    std::string hdr_str = "@HD	VN:1.0	SO:unsorted\n@SQ	SN:chr1	LN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;

    // Add 3 BAM_CDELs in the middle
    kstring_t str;
    std::string record = "16M_10M3D6M	0	chr1	674820	42	16M	*	0	0	ATTACATTCCATTCCA	~~~~~~~~~~~~~~~~";
    str.s = (char*) record.c_str();
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 3, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 6, BAM_CMATCH));
}


TEST(ChainTest, LiftCigar5) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    std::string hdr_str = "@HD	VN:1.0	SO:unsorted\n@SQ	SN:chr1	LN:248956422";
    sam_hdr_t* sam_hdr = sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;

    // Add 3 BAM_CSOFT_CLIPs in the beginning (secondary alignment)
    kstring_t str;
    std::string record = "16M_3S13M	256	chr1	687455	42	16M	*	0	0	*	*";
    str.s = (char*) record.c_str();
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CMATCH));
}


int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
