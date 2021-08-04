#include <iostream>
#include <unistd.h>
#include <string>
#include "leviosam.hpp"
#include "chain.hpp"
#include "gtest/gtest.h"

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
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
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
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
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
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
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
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
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
TEST(ChainMap, SimpleRankAndLift) {
    chain::ChainMap cmap ("small.chain", 0);

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


TEST(ChainMap, ParseChainLineCornerZero) {
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


TEST(ChainMap, LiftInReversedRegion) {
    chain::ChainMap cmap ("chr1_reversed_region.chain", 0);

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

TEST(ChainMap, LiftCigar) {
    chain::ChainMap cmap ("small.chain", 0);
    samFile* sam_fp = sam_open("chain_cigar.sam", "r");
    bam_hdr_t* sam_hdr = sam_hdr_read(sam_fp);
    bam1_t* aln = bam_init1();
    int err;
    size_t x;
    uint32_t* test_cigar;
    // Note: can use the helper function to print out CIGAR results
    // chain::debug_print_cigar(aln);

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
}

