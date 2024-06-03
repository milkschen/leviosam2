/*
 * leviosam_utils_test.cpp
 *
 * Test levioSAM utility functions
 *
 * Authors: Taher Mun, Nae-Chyun Chen, Ben Langmead
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "leviosam_utils.hpp"

#include <unistd.h>

#include <iostream>
#include <string>

#include "bed.hpp"
#include "gtest/gtest.h"
#include "htslib/sam.h"
#include "leviosam.hpp"

/* Bit_vector tests */

typedef sdsl::bit_vector::size_type size_type;

/**
 * Gets a fake sam header for chr1.
 *
 * @return The sam header
 */
bam_hdr_t* get_chr1_fake_hdr() {
    std::string hdr_str =
        "@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:248956422";
    return sam_hdr_parse(hdr_str.length(), &hdr_str[0]);
}

/**
 * Gets a fake bam record in chr1.
 * Uses chr1 to be compatible with get_chr1_fake_hdr(). Can generalize to other
 * chromosomes if needed.
 *
 * @param read The read string
 * @param sam_hdr The sam header
 * @param err The error code
 * @return The bam record
 */
bam1_t* get_chr1_fake_aln(const std::string read, sam_hdr_t* sam_hdr,
                          int& err) {
    bam1_t* aln = bam_init1();
    kstring_t kstr;
    kstr.s = (char*)read.c_str();
    kstr.l = read.length();
    kstr.m = kstr_get_m(kstr.l);

    err = sam_parse1(&kstr, sam_hdr, aln);
    return aln;
}

TEST(BitVector, SimpleCreateAndAccess) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    EXPECT_EQ(0, bv[0]);
    EXPECT_EQ(0, bv[500]);
    EXPECT_EQ(0, bv[10000 - 1]);
}

TEST(BitVector, WriteAndRead) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    EXPECT_EQ(1, bv[0]);
    EXPECT_EQ(0, bv[1]);
    EXPECT_EQ(0, bv[499]);
    EXPECT_EQ(1, bv[500]);
    EXPECT_EQ(0, bv[501]);
    EXPECT_EQ(0, bv[10000 - 2]);
    EXPECT_EQ(1, bv[10000 - 1]);
}

TEST(BitVector, RankSupport) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    sdsl::rank_support_v<> rs(&bv);
    EXPECT_EQ(0, rs.rank(0));
    EXPECT_EQ(1, rs.rank(1));
    EXPECT_EQ(1, rs.rank(2));
    EXPECT_EQ(1, rs.rank(10));
    EXPECT_EQ(1, rs.rank(500));
    EXPECT_EQ(2, rs.rank(501));
    EXPECT_EQ(2, rs.rank(601));
    EXPECT_EQ(2, rs.rank(10000 - 1));
    EXPECT_EQ(3, rs.rank(10000 - 0));
}

TEST(BitVector, SelectSupport) {
    sdsl::bit_vector bv(10000, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    sdsl::bit_vector::select_1_type ss(&bv);
    EXPECT_EQ(0, ss.select(1));
    EXPECT_EQ(500, ss.select(2));
    EXPECT_EQ(10000 - 1, ss.select(3));
}

TEST(CompressedVector, CreateAndRead) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    sdsl::sd_vector<> sv{bv};
    EXPECT_EQ(1, sv[0]);
    EXPECT_EQ(0, sv[1]);
    EXPECT_EQ(0, sv[499]);
    EXPECT_EQ(1, sv[500]);
    EXPECT_EQ(0, sv[501]);
    EXPECT_EQ(0, sv[10000 - 2]);
    EXPECT_EQ(1, sv[10000 - 1]);

    sdsl::rrr_vector<63> rv{bv};
    EXPECT_EQ(1, rv[0]);
    EXPECT_EQ(0, rv[1]);
    EXPECT_EQ(0, rv[499]);
    EXPECT_EQ(1, rv[500]);
    EXPECT_EQ(0, rv[501]);
    EXPECT_EQ(0, rv[10000 - 2]);
    EXPECT_EQ(1, rv[10000 - 1]);
}

TEST(CompressedVector, RankSupport) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    sdsl::sd_vector<> sv{bv};
    sdsl::sd_vector<>::rank_1_type rs(&sv);
    EXPECT_EQ(0, rs.rank(0));
    EXPECT_EQ(1, rs.rank(1));
    EXPECT_EQ(1, rs.rank(2));
    EXPECT_EQ(1, rs.rank(10));
    EXPECT_EQ(1, rs.rank(500));
    EXPECT_EQ(2, rs.rank(501));
    EXPECT_EQ(2, rs.rank(601));
    EXPECT_EQ(2, rs.rank(10000 - 1));
    EXPECT_EQ(3, rs.rank(10000 - 0));
}

TEST(CompressedVector, SelectSupport) {
    sdsl::bit_vector bv(10000, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000 - 1] = 1;
    sdsl::sd_vector<> sv{bv};
    sdsl::sd_vector<>::select_1_type ss{&sv};
    EXPECT_EQ(0, ss.select(1));
    EXPECT_EQ(500, ss.select(2));
    EXPECT_EQ(10000 - 1, ss.select(3));
}

/* Lift tests */

TEST(Lift, EmptyLift) {
    sdsl::bit_vector ibv;
    sdsl::bit_vector dbv;
    sdsl::bit_vector sbv;  // ins, del, snp
    ibv.resize(20);
    dbv.resize(20);
    sbv.resize(20);
    lift::Lift lift(ibv, dbv, sbv);
    EXPECT_EQ(lift.lift_pos(5), 5);
}

TEST(Lift, SimpleInsLift) {
    sdsl::bit_vector ibv;
    sdsl::bit_vector dbv;
    sdsl::bit_vector sbv;  // ins, del, snp
    ibv.resize(20);
    dbv.resize(20);
    sbv.resize(20);
    ibv[5] = 1;
    //     01234 56
    // I : 00000100000000000000
    // D : 00000000000000000000
    //     01234567
    lift::Lift lift(ibv, dbv, sbv);
    EXPECT_EQ(lift.lift_pos(7), 6);
}

TEST(Lift, SimpleDelLift) {
    sdsl::bit_vector ibv;
    sdsl::bit_vector dbv;
    sdsl::bit_vector sbv;  // ins, del, snp
    ibv.resize(20);
    dbv.resize(20);
    sbv.resize(20);
    dbv[5] = 1;
    //     01234567
    // I : 00000000000000000000
    // D : 00000100000000000000
    //     01234 56
    lift::Lift lift(ibv, dbv, sbv);
    EXPECT_EQ(lift.lift_pos(6), 7);
}

TEST(LiftMap, SimpleBamLift) {
    vcfFile* vcf_fp = bcf_open("simple_example.vcf", "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    lift::LiftMap lmap(vcf_fp, vcf_hdr, "", "0");
    samFile* sam_fp = sam_open("simple_example.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    err = sam_read1(sam_fp, sam_hdr, aln);
    x = lmap.lift_pos(sam_hdr->target_name[aln->core.tid], aln->core.pos);
    EXPECT_EQ(x, 5313299);
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    x = lmap.lift_pos(sam_hdr->target_name[aln->core.tid], aln->core.pos);
    EXPECT_EQ(x, 5315929);

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(LiftMap, SimpleBamCigarLift) {
    vcfFile* vcf_fp = bcf_open("simple_example.vcf", "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    lift::LiftMap lmap(vcf_fp, vcf_hdr, "", "0");
    samFile* sam_fp = sam_open("simple_example.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;
    err = sam_read1(sam_fp, sam_hdr, aln);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(11, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(2, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(9, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(15, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(1, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(5, BAM_CMATCH));

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(LiftMap, DelInIndelBamCigarLift) {
    vcfFile* vcf_fp = bcf_open("del_in_indel_example.vcf", "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    lift::LiftMap lmap(vcf_fp, vcf_hdr, "", "0");
    samFile* sam_fp = sam_open("del_in_indel_example.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;
    err = sam_read1(sam_fp, sam_hdr, aln);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(15, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(13, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(8, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(10, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(8, BAM_CMATCH));

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(LiftMap, SimpleSplicedBamCigarLift) {
    vcfFile* vcf_fp = bcf_open("spliced_example.vcf", "r");
    bcf_hdr_t* vcf_hdr = bcf_hdr_read(vcf_fp);
    lift::LiftMap lmap(vcf_fp, vcf_hdr, "", "0");
    samFile* sam_fp = sam_open("spliced_example.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;
    err = sam_read1(sam_fp, sam_hdr, aln);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(15, BAM_CREF_SKIP));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(13, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen(8, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen(10, BAM_CREF_SKIP));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen(8, BAM_CMATCH));

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(HtslibTest, CigarOpType) {
    auto typeM = bam_cigar_type(BAM_CMATCH);
    auto typeD = bam_cigar_type(BAM_CDEL);
    auto typeI = bam_cigar_type(BAM_CINS);
    auto typeS = bam_cigar_type(BAM_CSOFT_CLIP);
    EXPECT_EQ(typeM & 1, 1);
    EXPECT_EQ(typeM & 2, 2);
    EXPECT_EQ(typeM & 3, 3);
    EXPECT_EQ(bool(typeM & 2), true);

    EXPECT_EQ(typeD & 1, 0);
    EXPECT_EQ(typeD & 2, 2);
    EXPECT_EQ(typeD & 3, 2);
    EXPECT_EQ(typeD, 2);

    EXPECT_EQ(typeI & 1, 1);
    EXPECT_EQ(typeI & 2, 0);
    EXPECT_EQ(typeI & 3, 1);

    EXPECT_EQ(typeS & 1, 1);
    EXPECT_EQ(typeS & 2, 0);
    EXPECT_EQ(typeS & 3, 1);
}

TEST(UtilsTest, StrToSet) {
    std::set<std::string> s = LevioSamUtils::str_to_set("a,b", ",");
    EXPECT_EQ(s.size(), 2);
    ASSERT_NE(s.find("a"), s.end());
    ASSERT_NE(s.find("b"), s.end());
    EXPECT_EQ(s.find("c"), s.end());
}

TEST(UtilsTest, UpdateFlagUnmap) {
    samFile* sam_fp = sam_open("simple_example.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    LevioSamUtils::update_flag_unmap(aln, true, false);
    EXPECT_EQ(aln->core.qual, 0);

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    LevioSamUtils::update_flag_unmap(aln, true, true);
    EXPECT_EQ(aln->core.qual, 60);

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(UtilsTest, CheckSplitRule) {
    EXPECT_EQ(LevioSamUtils::check_split_rule("lifted"), true);
    EXPECT_EQ(LevioSamUtils::check_split_rule("mapq"), true);
    EXPECT_EQ(LevioSamUtils::check_split_rule("wrong_rule"), false);
}

TEST(UtilsTest, AddSplitRule) {
    LevioSamUtils::SplitRules split_rules;
    EXPECT_EQ(LevioSamUtils::add_split_rule(split_rules, "mapq:20"), true);
    EXPECT_EQ(split_rules[0].first, "mapq");
    EXPECT_FLOAT_EQ(split_rules[0].second, 20.);

    EXPECT_EQ(LevioSamUtils::add_split_rule(split_rules, "lifted"), true);
    EXPECT_EQ(split_rules[1].first, "lifted");
    EXPECT_FLOAT_EQ(split_rules[1].second, 1.);

    EXPECT_EQ(LevioSamUtils::add_split_rule(split_rules, "clipped_frac:0.66"),
              true);
    EXPECT_EQ(split_rules[2].first, "clipped_frac");
    EXPECT_FLOAT_EQ(split_rules[2].second, 0.66);

    EXPECT_EQ(LevioSamUtils::add_split_rule(split_rules, "none"), false);
    EXPECT_EQ(LevioSamUtils::add_split_rule(split_rules, "mapq:20:20"), false);
}

TEST(UltimaGenomicsTest, UpdateFlags) {
    samFile* sam_fp = sam_open("ultima_small.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    LevioSamUtils::update_ultima_genomics_tags(aln, true);

    // Tests the tp tag
    uint8_t* tp_ptr = bam_aux_get(aln, "tp");
    EXPECT_EQ(bam_auxB_len(tp_ptr), 10);
    std::vector<int8_t> exp_tp_array{-1, -1, 0, 1, -1, 1, 1, 1, 1, -1};
    for (uint32_t i = 0; i < 10; i++) {
        EXPECT_EQ(bam_auxB2i(tp_ptr, i), exp_tp_array[i]);
    }

    // Tests the t0 tag
    uint8_t* t0_ptr = bam_aux_get(aln, "t0");
    std::string t0_str = bam_aux2Z(t0_ptr);
    EXPECT_EQ(t0_str, "IIIII22222");

    sam_hdr_destroy(sam_hdr);
    sam_close(sam_fp);
}

TEST(ChainTest, GetMateQueryLenOnRef) {
    std::string record1 =
        "read1\t81\tchr1\t145334831\t33\t6S10M\t"
        "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:10M";
    int err;
    sam_hdr_t* sam_hdr = get_chr1_fake_hdr();
    bam1_t* aln = get_chr1_fake_aln(record1, sam_hdr, err);
    EXPECT_EQ(err, 0);
    EXPECT_EQ(LevioSamUtils::get_mate_query_len_on_ref(aln), 10);
    bam_destroy1(aln);

    std::string record2 =
        "read1\t81\tchr1\t145334831\t33\t6S10M\t"
        "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:6S5M1I4M";
    aln = get_chr1_fake_aln(record2, sam_hdr, err);
    EXPECT_EQ(err, 0);
    EXPECT_EQ(LevioSamUtils::get_mate_query_len_on_ref(aln), 9);
    bam_destroy1(aln);

    bam_hdr_destroy(sam_hdr);
}

TEST(WriteDeferredTest, WriteDeferredInit) {
    LevioSamUtils::WriteDeferred wd;
    LevioSamUtils::SplitRules split_rules;
    sam_hdr_t* hdr_orig;
    sam_hdr_t* hdr;
    BedUtils::Bed bed_defer_source, bed_defer_dest, bed_commit_source,
        bed_commit_dest;

    wd.init(split_rules, hdr_orig, hdr, bed_defer_source, bed_defer_dest,
            bed_commit_source, bed_commit_dest, 0);

    EXPECT_EQ(wd.get_write_deferred(), true);
}

TEST(WriteDeferredTest, CommitAlnDestMappingQuality) {
    LevioSamUtils::WriteDeferred wd;
    LevioSamUtils::SplitRules split_rules;
    LevioSamUtils::add_split_rule(split_rules, "mapq:30");

    sam_hdr_t* hdr_orig = NULL;  // not used in this test, so set to NULL.
    sam_hdr_t* hdr = get_chr1_fake_hdr();
    BedUtils::Bed bed_defer_source, bed_defer_dest, bed_commit_source,
        bed_commit_dest;
    wd.init(split_rules, hdr_orig, hdr, bed_defer_source, bed_defer_dest,
            bed_commit_source, bed_commit_dest, 0);

    int err;

    // Commits when MAPQ >= threshold
    std::string record_high_mapq =
        "read1\t81\tchr1\t145334831\t30\t16M\t"
        "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:6S5M1I4M";
    bam1_t* aln = get_chr1_fake_aln(record_high_mapq, hdr, err);

    EXPECT_EQ(err, 0);
    EXPECT_EQ(wd.commit_aln_dest(aln), true);
    bam_destroy1(aln);

    // defers when MAPQ < threshold
    std::string record_low_mapq =
        "read1\t81\tchr1\t145334831\t29\t16M\t"
        "=\t1245932\t0\tATTACATTCCATTCCA\t~~~~~~~~~~~~~~~~\tMC:Z:6S5M1I4M";
    aln = get_chr1_fake_aln(record_low_mapq, hdr, err);
    EXPECT_EQ(err, 0);
    EXPECT_EQ(wd.commit_aln_dest(aln), false);
    bam_destroy1(aln);

    bam_hdr_destroy(hdr);
}

TEST(FastqRecordTest, Init) {
    LevioSamUtils::FastqRecord fr("ATCG", "CCC!");
    EXPECT_EQ(fr.seq_str, "ATCG");
    EXPECT_EQ(fr.qual_str, "CCC!");
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
