/*
 * reconcile_test.cpp
 *
 * Test the reconcile module of levioSAM2
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "reconcile.hpp"

#include <unistd.h>

#include <iostream>
#include <vector>

#include "gtest/gtest.h"
#include "leviosam_utils.hpp"

/* Reconcile tests */
TEST(ReconcileTest, SelectBestAlnOneBestAtFirst) {
    std::vector<bool> pi = {true, true, false};
    std::vector<int> sc = {0, -3, -20};
    std::vector<int> mq = {40, 30, 2};
    int ntb = 0;

    int rank = select_best_aln(pi, sc, mq, ntb);

    EXPECT_EQ(rank, 0);
    EXPECT_EQ(ntb, 1);
}

TEST(ReconcileTest, SelectBestAlnOneBestAtSecond) {
    std::vector<bool> pi = {true, true, false};
    std::vector<int> sc = {0, 3, -20};
    std::vector<int> mq = {40, 30, 2};
    int ntb = 0;

    int rank = select_best_aln(pi, sc, mq, ntb);

    EXPECT_EQ(rank, 1);
    EXPECT_EQ(ntb, 1);
}

TEST(ReconcileTest, SelectBestAlnFourBest) {
    std::vector<bool> pi = {true, true, true, true, false};
    std::vector<int> sc = {0, 0, 0, 0, -20};
    std::vector<int> mq = {40, 40, 40, 40, 2};
    int ntb = 0;

    int rank = select_best_aln(pi, sc, mq, ntb);

    EXPECT_EQ(rank, 2);
    EXPECT_EQ(ntb, 4);
}

