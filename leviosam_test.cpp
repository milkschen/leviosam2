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
    std::vector<std::string> chroms_not_found;
    std::mutex mutex_vec;
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
    std::vector<std::string> chroms_not_found;
    std::mutex mutex_vec;
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
    std::vector<std::string> chroms_not_found;
    std::mutex mutex_vec;
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
    std::vector<std::string> chroms_not_found;
    std::mutex mutex_vec;
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

    std::vector<std::string> uc;
    std::mutex mutex;

    int pos = 674047;
    std::string contig = "chr1";
    int rank = cmap.get_start_rank(contig, pos);
    // std::cerr << "rank(" << pos << ")=" << rank << "\n";
    EXPECT_EQ(rank, 0);
    EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
    EXPECT_EQ(cmap.lift_pos(contig, pos), 100272);
    

    pos = 207130;
    contig = "chr2";
    rank = cmap.get_start_rank(contig, pos);
    // std::cerr << "rank(" << pos << ")=" << rank << "\n";
    EXPECT_EQ(rank, 1);
    EXPECT_EQ(cmap.lift_contig(contig, pos), contig);
    EXPECT_EQ(cmap.lift_pos(contig, pos), 206846+121+8);
}
