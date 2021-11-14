#include <iostream>
#include <unistd.h>
#include <string>
#include "chain.hpp"
#include "leviosam.hpp"
#include "leviosam_utils.hpp"
#include "htslib/sam.h"
#include "gtest/gtest.h"
#include "bed.hpp"

//
// bit_vector tests
//

typedef sdsl::bit_vector::size_type size_type;

TEST(BitVector, SimpleCreateAndAccess) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    EXPECT_EQ(0, bv[0]);
    EXPECT_EQ(0, bv[500]);
    EXPECT_EQ(0, bv[10000-1]);
}

TEST(BitVector, WriteAndRead) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    EXPECT_EQ(1, bv[0]);
    EXPECT_EQ(0, bv[1]);
    EXPECT_EQ(0, bv[499]);
    EXPECT_EQ(1, bv[500]);
    EXPECT_EQ(0, bv[501]);
    EXPECT_EQ(0, bv[10000-2]);
    EXPECT_EQ(1, bv[10000-1]);
}

TEST(BitVector, RankSupport) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    sdsl::rank_support_v<> rs(&bv);
    EXPECT_EQ(0, rs.rank(0));
    EXPECT_EQ(1, rs.rank(1));
    EXPECT_EQ(1, rs.rank(2));
    EXPECT_EQ(1, rs.rank(10));
    EXPECT_EQ(1, rs.rank(500));
    EXPECT_EQ(2, rs.rank(501));
    EXPECT_EQ(2, rs.rank(601));
    EXPECT_EQ(2, rs.rank(10000-1));
    EXPECT_EQ(3, rs.rank(10000-0));
}

TEST(BitVector, SelectSupport) {
    sdsl::bit_vector bv(10000, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    sdsl::bit_vector::select_1_type ss(&bv);
    EXPECT_EQ(0, ss.select(1));
    EXPECT_EQ(500, ss.select(2));
    EXPECT_EQ(10000-1, ss.select(3));
}

TEST(CompressedVector, CreateAndRead) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    sdsl::sd_vector<> sv{bv};
    EXPECT_EQ(1, sv[0]);
    EXPECT_EQ(0, sv[1]);
    EXPECT_EQ(0, sv[499]);
    EXPECT_EQ(1, sv[500]);
    EXPECT_EQ(0, sv[501]);
    EXPECT_EQ(0, sv[10000-2]);
    EXPECT_EQ(1, sv[10000-1]);

    sdsl::rrr_vector<63> rv{bv};
    EXPECT_EQ(1, rv[0]);
    EXPECT_EQ(0, rv[1]);
    EXPECT_EQ(0, rv[499]);
    EXPECT_EQ(1, rv[500]);
    EXPECT_EQ(0, rv[501]);
    EXPECT_EQ(0, rv[10000-2]);
    EXPECT_EQ(1, rv[10000-1]);
}

TEST(CompressedVector, RankSupport) {
    size_type n = 10000;
    sdsl::bit_vector bv(n, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    sdsl::sd_vector<> sv{bv};
    sdsl::sd_vector<>::rank_1_type rs(&sv);
    EXPECT_EQ(0, rs.rank(0));
    EXPECT_EQ(1, rs.rank(1));
    EXPECT_EQ(1, rs.rank(2));
    EXPECT_EQ(1, rs.rank(10));
    EXPECT_EQ(1, rs.rank(500));
    EXPECT_EQ(2, rs.rank(501));
    EXPECT_EQ(2, rs.rank(601));
    EXPECT_EQ(2, rs.rank(10000-1));
    EXPECT_EQ(3, rs.rank(10000-0));
}

TEST(CompressedVector, SelectSupport) {
    sdsl::bit_vector bv(10000, 0);
    bv[0] = 1;
    bv[500] = 1;
    bv[10000-1] = 1;
    sdsl::sd_vector<> sv{bv};
    sdsl::sd_vector<>::select_1_type ss{&sv};
    EXPECT_EQ(0, ss.select(1));
    EXPECT_EQ(500, ss.select(2));
    EXPECT_EQ(10000-1, ss.select(3));
}

/* Lift tests */

TEST(Lift, EmptyLift) {
    sdsl::bit_vector ibv;
    sdsl::bit_vector dbv;
    sdsl::bit_vector sbv; // ins, del, snp
    ibv.resize(20);
    dbv.resize(20);
    sbv.resize(20);
    lift::Lift lift(ibv, dbv, sbv);
    EXPECT_EQ(lift.lift_pos(5), 5);
}

TEST(Lift, SimpleInsLift) {
    sdsl::bit_vector ibv;
    sdsl::bit_vector dbv;
    sdsl::bit_vector sbv; // ins, del, snp
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
    sdsl::bit_vector sbv; // ins, del, snp
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 11, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 2, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 9, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 15, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 1, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 5, BAM_CMATCH));
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 15, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 13, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 8, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 10, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 8, BAM_CMATCH));
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
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 15, BAM_CREF_SKIP));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 13, BAM_CMATCH));

    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    lmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 8, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 10, BAM_CREF_SKIP));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 8, BAM_CMATCH));
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
    EXPECT_EQ(cmap.lift_pos(contig, pos), 100272);

    pos = 207130;
    contig = "chr2";
    rank = cmap.get_start_rank(contig, pos);
    EXPECT_EQ(rank, 2);
    EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
    EXPECT_EQ(cmap.lift_pos(contig, pos), 206846+121+8);
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
        EXPECT_EQ(cmap.lift_pos(contig, pos), gold_pos_array[i]);
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
    LevioSamUtils::debug_print_cigar(new_cigar.data(), new_cigar.size());
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


TEST(ChainTest, LiftCigar) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap cmap ("small.chain", 0, 0, lm);
    samFile* sam_fp = sam_open("chain_cigar.sam", "r");
    sam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;
    // Note: can use the helper function to print out CIGAR results
    // chain::debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);

    // CIGAR should not change
    err = sam_read1(sam_fp, sam_hdr, aln);
    cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 7, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 13, BAM_CMATCH));

    // Add 3 BAM_CSOFT_CLIPs in the beginning
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CMATCH));
    
    // Add 2 BAM_CINSs in the middle
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 2, BAM_CINS));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 4, BAM_CMATCH));
    
    // Add 3 BAM_CDELs in the middle
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 3);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 10, BAM_CMATCH));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 3, BAM_CDEL));
    EXPECT_EQ(test_cigar[2], bam_cigar_gen( 6, BAM_CMATCH));

    // Add 3 BAM_CSOFT_CLIPs in the beginning (secondary alignment)
    err = sam_read1(sam_fp, sam_hdr, aln);
    EXPECT_EQ(err, 0);
    cmap.lift_cigar(sam_hdr->target_name[aln->core.tid], aln);
    test_cigar = bam_get_cigar(aln);
    EXPECT_EQ(aln->core.n_cigar, 2);
    EXPECT_EQ(test_cigar[0], bam_cigar_gen( 3, BAM_CSOFT_CLIP));
    EXPECT_EQ(test_cigar[1], bam_cigar_gen( 13, BAM_CMATCH));
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


TEST(BedTest, AddBedRecord) {
    auto b = BedUtils::Bed();
    b.add_interval("chr1	972	1052");
    b.add_interval("chr1	1053	1085");
    b.add_interval("chr1	1150	1198");
    b.add_interval("chr1	1199	1262");
    b.add_interval("chr2	340	498");
    b.add_interval("chr2	599	1000");
    int idx_contig_cnt = b.index();
    int contig_cnt = 0;
    // Index the interval trees for each contig
    for (auto& itr: b.get_intervals()) {
        contig_cnt += 1;
    }
    EXPECT_EQ(contig_cnt, 2);
    EXPECT_EQ(contig_cnt, idx_contig_cnt);
    EXPECT_EQ(b.intersect("chr1", 1000), 1);
    EXPECT_EQ(b.intersect("chr1", 1052), 0);
    EXPECT_EQ(b.intersect("chr2", 339), 0);
    EXPECT_EQ(b.intersect("chr2", 340), 1);
    EXPECT_EQ(b.intersect("chr2", 497), 1);
    EXPECT_EQ(b.intersect("chr2", 498), 0);
    EXPECT_EQ(b.intersect("chr2", 340), 1);
    EXPECT_EQ(b.intersect("chr2", 340), 1);
    EXPECT_EQ(b.intersect("chr2", 598), 0);
    EXPECT_EQ(b.intersect("chr2", 599), 1);
    EXPECT_EQ(b.intersect("chr2", 600), 1);
    EXPECT_EQ(b.intersect("chr1", 1500), 0);
    EXPECT_EQ(b.intersect("chr15", 1000), 0);
}


TEST(UtilsTest, StrToSet) {
    std::set<std::string> s = LevioSamUtils::str_to_set("a,b", ",");
    EXPECT_EQ(s.size(), 2);
    ASSERT_NE(s.find("a"), s.end());
    ASSERT_NE(s.find("b"), s.end());
    EXPECT_EQ(s.find("c"), s.end());
}

int main(int argc, char **argv){
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
