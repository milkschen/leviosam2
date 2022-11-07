/*
 * bed_test.cpp
 *
 * Test bed related functions for levioSAM2
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "bed.hpp"

#include <unistd.h>

#include <iostream>

#include "gtest/gtest.h"
#include "leviosam_utils.hpp"

/* Bed tests */
TEST(BedTest, Simple) {
    std::string bed_fn("test_small.bed");
    BedUtils::Bed bed(bed_fn);
    EXPECT_EQ(bed.get_fn(), bed_fn);
    EXPECT_EQ(bed.index(), 3);
    EXPECT_EQ(bed.intersect("chr1", 450, 499), false);
    EXPECT_EQ(bed.intersect("chr1", 450, 501), true);
}

TEST(BedTest, AddInterval) {
    std::string bed_fn("test_small.bed");
    BedUtils::Bed bed(bed_fn);
    EXPECT_EQ(bed.intersect("chr1", 55), false);
    EXPECT_EQ(bed.index(), 3);
    EXPECT_EQ(bed.add_interval("chr1	50	60"), true);
    // chr1 already exists, so cnt is still 3
    EXPECT_EQ(bed.index(), 3);
    EXPECT_EQ(bed.intersect("chr1", 55), true);
    EXPECT_EQ(bed.intersect("chr1", 65), false);
    EXPECT_EQ(bed.intersect("chr1", 55, 58), true);
    EXPECT_EQ(bed.intersect("chr1", 45, 80), true);
    EXPECT_EQ(bed.intersect("chr1", 45, 55), true);
}

TEST(BedTest, IntersectWithFrac) {
    std::string bed_fn("test_small.bed");
    BedUtils::Bed bed(bed_fn);
    EXPECT_EQ(bed.intersect("chr1", 450, 650, 1.0), false);
    EXPECT_EQ(bed.intersect("chr1", 550, 580, 1.0), true);
    EXPECT_EQ(bed.intersect("chr1", 550, 580, 0.99), true);
    EXPECT_EQ(bed.intersect("chr1", 550, 700, 0.3), true);
    EXPECT_EQ(bed.intersect("chr1", 550, 700, 0.35), false);
    EXPECT_EQ(bed.intersect("chr1", 550, 700, 30), true);
    EXPECT_EQ(bed.intersect("chr1", 550, 580, 0), true);
    EXPECT_EQ(bed.intersect("chr1", 550, 580, -1), true);
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

