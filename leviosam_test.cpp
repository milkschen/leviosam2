#include "leviosam.hpp"
#include "gtest/gtest.h"

//
// bit_veector tests
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
