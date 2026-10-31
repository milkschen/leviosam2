#include <gtest/gtest.h>
#include <htslib/kstring.h>

#include <memory>
#include <cstring>
#include <stdexcept>
#include <string>
#include <utility>

#include "mate_metadata.hpp"

namespace {

using Bam = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

class MateMetadataTest : public ::testing::Test {
   protected:
    sam_hdr_t* header = nullptr;
    void SetUp() override {
        const char* text = "@HD\tVN:1.6\n@SQ\tSN:chr\tLN:1000000\n"
                           "@RG\tID:a\tSM:sample\tPL:ILLUMINA\n"
                           "@RG\tID:b\tSM:sample\tPL:ILLUMINA\n";
        header = sam_hdr_parse(std::strlen(text), text);
        ASSERT_NE(header, nullptr);
    }
    void TearDown() override { sam_hdr_destroy(header); }
    Bam parse(const std::string& name, int flag, int pos = 101,
               int mate_pos = 201, const std::string& group = "a") {
        Bam bam(bam_init1(), &bam_destroy1);
        std::string text = name + "\t" + std::to_string(flag) + "\tchr\t" +
            std::to_string(pos) + "\t60\t10M\t=\t" + std::to_string(mate_pos) +
            "\t0\tAAAAAAAAAA\tIIIIIIIIII\tRG:Z:" + group + "\tLO:Z:L_UM";
        kstring_t line = KS_INITIALIZE;
        kputs(text.c_str(), &line);
        int status = sam_parse1(&line, header, bam.get());
        free(line.s);
        EXPECT_GE(status, 0);
        return bam;
    }
};

TEST_F(MateMetadataTest, DiskIndexGrowsAndReplaysAllRecords) {
    // Start with four slots, forcing repeated growth and colliding probe
    // sequences without needing tens of thousands of fixture alignments.
    LevioSamUtils::MateMetadataSpool spool(4);
    const int count = 2048;
    for (int i = 0; i < count; ++i) {
        Bam source = parse("read" + std::to_string(i), 65);
        Bam lifted(bam_dup1(source.get()), &bam_destroy1);
        lifted->core.pos = 500 + i;
        spool.append(source.get(), lifted.get(), header, header);
    }
    for (int i = count - 1; i >= 0; --i) {
        Bam mate = parse("read" + std::to_string(i), 129);
        LevioSamUtils::MateOutcome result;
        ASSERT_TRUE(spool.find_mate(mate.get(), result));
        EXPECT_EQ(result.pos, 500 + i);
        EXPECT_EQ(result.end, 510 + i);
    }
    spool.begin_replay(header, header);
    Bam source(bam_init1(), &bam_destroy1), lifted(bam_init1(), &bam_destroy1);
    for (int i = 0; i < count; ++i) {
        ASSERT_TRUE(spool.next(source.get(), lifted.get(), header, header));
        EXPECT_EQ(std::string(bam_get_qname(source.get())), "read" + std::to_string(i));
        EXPECT_EQ(source->core.pos, 100);
        EXPECT_EQ(lifted->core.pos, 500 + i);
    }
    EXPECT_FALSE(spool.next(source.get(), lifted.get(), header, header));
}

TEST_F(MateMetadataTest, ReadGroupsSeparateIdenticalNamesAndDuplicatesFail) {
    LevioSamUtils::MateMetadataSpool spool(4);
    Bam source_a = parse("same", 65, 101, 201, "a");
    Bam source_b = parse("same", 65, 301, 401, "b");
    spool.append(source_a.get(), source_a.get(), header, header);
    spool.append(source_b.get(), source_b.get(), header, header);
    for (const auto& group : {"a", "b"}) {
        Bam opposite = parse("same", 129, 101, 201, group);
        LevioSamUtils::MateOutcome result;
        ASSERT_TRUE(spool.find_mate(opposite.get(), result));
        EXPECT_EQ(result.pos, std::string(group) == "a" ? 100 : 300);
    }
    EXPECT_THROW(spool.append(source_a.get(), source_a.get(), header, header),
                 std::runtime_error);
    Bam absent = parse("absent", 129);
    LevioSamUtils::MateOutcome result;
    EXPECT_FALSE(spool.find_mate(absent.get(), result));
}

TEST_F(MateMetadataTest, PreservesProperLibraryOrientationAndUsesActualSpan) {
    // Preserve FR, RF, FF and RR library conventions under translation.
    for (const auto& flags : {std::make_pair(99, 147), std::make_pair(83, 163),
                              std::make_pair(67, 131), std::make_pair(115, 179)}) {
        Bam source = parse("pair", flags.first, 101, 201);
        Bam mate_source = parse("pair", flags.second, 201, 101);
        Bam lifted = parse("pair", flags.first, 301, 401);
        Bam mate = parse("pair", flags.second, 401, 301);
        auto outcome = LevioSamUtils::mate_outcome(mate.get(), mate_source.get());
        LevioSamUtils::repair_mate_metadata(lifted.get(), source.get(), &outcome);
        EXPECT_TRUE(lifted->core.flag & BAM_FPROPER_PAIR) << flags.first;
        EXPECT_EQ(lifted->core.mpos, 400);
        EXPECT_EQ(lifted->core.isize, 110);
        EXPECT_EQ(bool(lifted->core.flag & BAM_FMREVERSE),
                   bool(mate->core.flag & BAM_FREVERSE));
    }
}

TEST_F(MateMetadataTest, CorrectsFlagsAndUnmappedPlacementFromActualMate) {
    Bam source = parse("pair", 99);
    Bam mate_source = parse("pair", 147, 201, 101);
    Bam lifted = parse("pair", 99, 301, 401);
    Bam mate = parse("pair", 133, 301, 301);
    auto outcome = LevioSamUtils::mate_outcome(mate.get(), mate_source.get());
    LevioSamUtils::repair_mate_metadata(lifted.get(), source.get(), &outcome);
    EXPECT_TRUE(lifted->core.flag & BAM_FMUNMAP);
    EXPECT_FALSE(lifted->core.flag & BAM_FMREVERSE);
    EXPECT_FALSE(lifted->core.flag & BAM_FPROPER_PAIR);
    EXPECT_EQ(lifted->core.mpos, lifted->core.pos);
    EXPECT_EQ(lifted->core.isize, 0);
    EXPECT_STREQ(bam_aux2Z(bam_aux_get(lifted.get(), "LO")), "L_UM");

    // A missing mate is unknown, not falsely marked unmapped.
    LevioSamUtils::repair_mate_metadata(lifted.get(), source.get(), nullptr);
    EXPECT_FALSE(lifted->core.flag & BAM_FMUNMAP);
    EXPECT_EQ(lifted->core.mtid, -1);
    EXPECT_EQ(lifted->core.mpos, -1);
    EXPECT_EQ(lifted->core.isize, 0);
    EXPECT_STREQ(bam_aux2Z(bam_aux_get(lifted.get(), "LO")), "L_L");
}

TEST_F(MateMetadataTest, SupplementaryProperStatusFollowsActualPrimaryPair) {
    Bam source = parse("pair", 2195, 501, 101);
    Bam lifted(bam_dup1(source.get()), &bam_destroy1);
    source->core.tid = 1;
    lifted->core.tid = 1;  // supplementary alignment on another chromosome
    Bam own_source = parse("pair", 147, 201, 101);
    Bam own_lifted = parse("pair", 147, 401, 301);
    Bam mate_source = parse("pair", 99, 101, 201);
    Bam mate_lifted = parse("pair", 99, 301, 401);
    auto own = LevioSamUtils::mate_outcome(own_lifted.get(), own_source.get());
    auto mate = LevioSamUtils::mate_outcome(mate_lifted.get(), mate_source.get());
    LevioSamUtils::repair_mate_metadata(lifted.get(), source.get(), &mate, &own);
    EXPECT_TRUE(lifted->core.flag & BAM_FPROPER_PAIR);
    EXPECT_EQ(lifted->core.mtid, 0);
    EXPECT_EQ(lifted->core.mpos, 300);
    EXPECT_EQ(lifted->core.isize, 0);
    // If the actual primary pair ceases to be proper, inherited status clears.
    mate.flag |= BAM_FUNMAP;
    LevioSamUtils::repair_mate_metadata(lifted.get(), source.get(), &mate, &own);
    EXPECT_FALSE(lifted->core.flag & BAM_FPROPER_PAIR);
}

}  // namespace
