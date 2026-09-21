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

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "gtest/gtest.h"
#include "leviosam_utils.hpp"

namespace {

using BamPtr = std::unique_ptr<bam1_t, decltype(&bam_destroy1)>;

std::string write_temporary_file(const std::string &contents) {
    char path[] = "/tmp/leviosam2-chain-matrix-XXXXXX";
    int fd = mkstemp(path);
    if (fd < 0) return "";
    close(fd);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out << contents;
    out.close();
    return path;
}

BamPtr parse_alignment(sam_hdr_t *header, const std::string &contig,
                       int32_t position, const std::string &cigar) {
    std::string record = "matrix\t0\t" + contig + "\t" +
                         std::to_string(position + 1) + "\t60\t" + cigar +
                         "\t*\t0\t0\t*\t*";
    kstring_t line = KS_INITIALIZE;
    kputs(record.c_str(), &line);
    BamPtr alignment(bam_init1(), bam_destroy1);
    if (sam_parse1(&line, header, alignment.get()) < 0) alignment.reset();
    free(line.s);
    return alignment;
}

std::string cigar_string(const bam1_t *alignment) {
    std::string result;
    const uint32_t *cigar = bam_get_cigar(alignment);
    for (uint32_t i = 0; i < alignment->core.n_cigar; ++i) {
        result += std::to_string(bam_cigar_oplen(cigar[i]));
        result += bam_cigar_opchr(cigar[i]);
    }
    return result;
}

struct CigarResult {
    int status;
    std::string cigar;
    hts_pos_t query_length;

    bool operator==(const CigarResult &other) const {
        return std::tie(status, cigar, query_length) ==
               std::tie(other.status, other.cigar, other.query_length);
    }
};

class SyntheticChainMatrixTest : public testing::Test {
   protected:
    void SetUp() override {
        static const std::string chains =
            "chain 1 src 100 + 0 100 dst 120 + 10 111 1\n"
            "20 0 3\n"
            "20 2 0\n"
            "20 0 0\n"
            "38\n"
            "\n"
            "chain 1 revsrc 100 + 0 100 revdst 120 - 10 108 2\n"
            "40 2 0\n"
            "58\n";
        chain_path = write_temporary_file(chains);
        ASSERT_FALSE(chain_path.empty());

        chain::LengthMap lengths{{"dst", 120}, {"revdst", 120}};
        parsed.reset(new chain::ChainMap(chain_path, 0, 10, lengths));

        index_path = write_temporary_file("");
        ASSERT_FALSE(index_path.empty());
        {
            std::ofstream out(index_path,
                              std::ios::binary | std::ios::trunc);
            ASSERT_TRUE(out.good());
            parsed->serialize(out);
        }
        std::ifstream in(index_path, std::ios::binary);
        ASSERT_TRUE(in.good());
        restored.reset(new chain::ChainMap(in, 0, 10));

        const std::string header_text =
            "@HD\tVN:1.6\tSO:unsorted\n"
            "@SQ\tSN:src\tLN:100\n"
            "@SQ\tSN:revsrc\tLN:100\n";
        header = sam_hdr_parse(header_text.size(), header_text.c_str());
        ASSERT_NE(header, nullptr);
    }

    void TearDown() override {
        sam_hdr_destroy(header);
        unlink(chain_path.c_str());
        unlink(index_path.c_str());
    }

    CigarResult lift_cigar(chain::ChainMap &chain_map,
                           const std::string &contig, int32_t position,
                           const std::string &cigar) {
        BamPtr alignment = parse_alignment(header, contig, position, cigar);
        EXPECT_NE(alignment, nullptr);
        if (!alignment) return {-1, "", 0};
        hts_pos_t query_length = bam_cigar2qlen(
            alignment->core.n_cigar, bam_get_cigar(alignment.get()));
        int status = chain_map.lift_cigar(contig, alignment.get());
        hts_pos_t lifted_query_length = bam_cigar2qlen(
            alignment->core.n_cigar, bam_get_cigar(alignment.get()));
        EXPECT_EQ(lifted_query_length, query_length);
        return {status, cigar_string(alignment.get()), lifted_query_length};
    }

    void expect_cigar(const std::string &contig, int32_t position,
                      const std::string &input, const std::string &expected) {
        CigarResult parsed_result =
            lift_cigar(*parsed, contig, position, input);
        CigarResult restored_result =
            lift_cigar(*restored, contig, position, input);
        EXPECT_EQ(parsed_result, restored_result);
        EXPECT_EQ(parsed_result.status, 0);
        EXPECT_EQ(parsed_result.cigar, expected);
    }

    std::string chain_path;
    std::string index_path;
    std::unique_ptr<chain::ChainMap> parsed;
    std::unique_ptr<chain::ChainMap> restored;
    sam_hdr_t *header = nullptr;
};

}  // namespace

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

TEST(ChainTest, SerializationRoundTrip) {
    std::vector<std::pair<std::string, int32_t>> lm;
    lm.push_back(std::make_pair("chr1", 248387328));
    chain::ChainMap original("small.chain", 0, 0, lm);

    char path[] = "/tmp/leviosam2-chain-roundtrip-XXXXXX";
    int fd = mkstemp(path);
    ASSERT_GE(fd, 0);
    close(fd);
    {
        std::ofstream out(path, std::ios::binary | std::ios::trunc);
        ASSERT_TRUE(out.good());
        original.serialize(out);
    }
    std::ifstream in(path, std::ios::binary);
    ASSERT_TRUE(in.good());
    chain::ChainMap restored(in, 0, 0);
    unlink(path);

    EXPECT_EQ(restored.length_map, original.length_map);
    EXPECT_EQ(restored.lift_contig("chr1", 674047), "chr1");
    EXPECT_EQ(restored.lift_pos("chr1", 674047, 0, true), 100272);
}

TEST_F(SyntheticChainMatrixTest, PositionBoundariesMatchAfterSerialization) {
    struct PositionCase {
        std::string contig;
        hts_pos_t source;
        hts_pos_t expected;
    };
    const std::vector<PositionCase> cases{
        {"src", 0, 10},  {"src", 19, 29}, {"src", 20, 33},
        {"src", 39, 52}, {"src", 40, 53}, {"src", 41, 53},
        {"src", 42, 53}, {"src", 99, 110}};

    for (const auto &test_case : cases) {
        SCOPED_TRACE(test_case.contig + ":" +
                     std::to_string(test_case.source));
        hts_pos_t parsed_position =
            parsed->lift_pos(test_case.contig, test_case.source, 2, true);
        hts_pos_t restored_position =
            restored->lift_pos(test_case.contig, test_case.source, 2, true);
        EXPECT_EQ(parsed_position, restored_position);
        EXPECT_EQ(parsed_position, test_case.expected);
        EXPECT_EQ(parsed->lift_contig(test_case.contig, test_case.source),
                  restored->lift_contig(test_case.contig, test_case.source));
    }

    EXPECT_EQ(parsed->lift_pos("missing", 10, 0, true), -1);
    EXPECT_EQ(restored->lift_pos("missing", 10, 0, true), -1);
    EXPECT_EQ(parsed->lift_contig("missing", 10), "*");
    EXPECT_EQ(restored->lift_contig("missing", 10), "*");
    EXPECT_EQ(parsed->lift_pos("src", 100, 0, true), -1);
    EXPECT_EQ(restored->lift_pos("src", 100, 0, true), -1);
    EXPECT_EQ(parsed->lift_contig("src", 100), "*");
    EXPECT_EQ(restored->lift_contig("src", 100), "*");
}

TEST_F(SyntheticChainMatrixTest, CigarBoundaryAndGapMatrix) {
    // Exact interval ends and starts do not acquire a chain-gap operation.
    expect_cigar("src", 15, "5M", "5M");
    expect_cigar("src", 20, "5M", "5M");

    // A target-only gap becomes a deletion; a source-only gap becomes an
    // insertion while preserving query length.
    expect_cigar("src", 15, "10M", "5M3D5M");
    expect_cigar("src", 35, "12M", "5M2I5M");
    expect_cigar("src", 15, "55M", "5M3D20M2I28M");
    expect_cigar("revsrc", 35, "12M", "5M2I5M");

    // Reads beginning or ending in a small internal source gap are clipped to
    // the neighboring mapped interval.
    expect_cigar("src", 41, "6M", "1S5M");
    expect_cigar("src", 35, "6M", "5M1S");
    expect_cigar("revsrc", 41, "6M", "5M1S");
    expect_cigar("revsrc", 35, "6M", "1S5M");
}

TEST_F(SyntheticChainMatrixTest, PreservesAndNormalizesCigarOperations) {
    expect_cigar("src", 5, "2S3=1X2I4M2D3N1H",
                 "2S4M2I4M2D3N1H");
    expect_cigar("revsrc", 20, "2S3=1X2I4M2D3N1H",
                 "1H3N2D4M2I4M2S");
}

TEST_F(SyntheticChainMatrixTest, PlacesGapsByReferenceConsumption) {
    // D and N consume the source reference but not the query. The chain
    // breakpoint must therefore be placed after them rather than after the
    // same number of query bases.
    expect_cigar("src", 10, "5M5D10M", "5M8D10M");
    expect_cigar("src", 10, "5M5N10M", "5M5N3D10M");

    // I consumes query but not source reference, so it must not advance the
    // chain breakpoint.
    expect_cigar("src", 10, "5M5I10M", "5M5I5M3D5M");
    expect_cigar("revsrc", 35, "3M2I9M", "5M2I2M2I3M");

    // A source-side chain deletion covered by a read deletion cancels that
    // deletion instead of manufacturing query-consuming inserted bases.
    expect_cigar("src", 30, "10M2D10M", "20M");
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

TEST(ChainTest, IntervalMapSanityCheckEmpty) {
    chain::ChainMap cmap;
    // Empty interval map should pass sanity check
    EXPECT_EQ(cmap.validate_intervals(), true);
}

TEST(ChainTest, IntervalMapSanityCheckSingleInterval) {
    chain::ChainMap cmap;
    // Single interval should pass sanity check
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    EXPECT_EQ(cmap.validate_intervals(), true);
}

TEST(ChainTest, IntervalMapSanityCheckNoOverlaps) {
    chain::ChainMap cmap;
    // Add non-overlapping intervals
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 250, 350, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 400, 500, 50, true));
    
    EXPECT_EQ(cmap.validate_intervals(), true);
}

TEST(ChainTest, IntervalMapSanityCheckWithOverlaps) {
    chain::ChainMap cmap;
    // Add overlapping intervals - first interval ends at 200, second starts at 150
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 150, 250, 50, true));
    
    EXPECT_EQ(cmap.validate_intervals(), false);
}

TEST(ChainTest, IntervalMapSanityCheckAdjacentIntervals) {
    chain::ChainMap cmap;
    // Add adjacent intervals (no gap, no overlap) - should pass
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 200, 300, 50, true));
    
    EXPECT_EQ(cmap.validate_intervals(), true);
}

TEST(ChainTest, IntervalMapSanityCheckMultipleContigs) {
    chain::ChainMap cmap;
    // Add intervals to multiple contigs - some with overlaps, some without
    // chr1: no overlaps (should pass)
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 250, 350, 50, true));
    
    // chr2: has overlaps (should fail)
    cmap.add_interval("chr2", chain::Interval("chr2_dest", 100, 200, 50, true));
    cmap.add_interval("chr2", chain::Interval("chr2_dest", 150, 250, 50, true));
    
    // Should fail because chr2 has overlaps
    EXPECT_EQ(cmap.validate_intervals(), false);
}

TEST(ChainTest, IntervalMapSanityCheckAllContigsValid) {
    chain::ChainMap cmap;
    // Add intervals to multiple contigs - all valid (no overlaps)
    // chr1: no overlaps
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 250, 350, 50, true));
    
    // chr2: no overlaps
    cmap.add_interval("chr2", chain::Interval("chr2_dest", 100, 200, 50, true));
    cmap.add_interval("chr2", chain::Interval("chr2_dest", 300, 400, 50, true));
    
    // chr3: single interval
    cmap.add_interval("chr3", chain::Interval("chr3_dest", 500, 600, 50, true));
    
    // Should pass because all contigs have valid intervals
    EXPECT_EQ(cmap.validate_intervals(), true);
}

TEST(ChainTest, IntervalMapSanityCheckExactOverlap) {
    chain::ChainMap cmap;
    // Add intervals with exact overlap (same start/end positions)
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 200, 50, true));
    
    EXPECT_EQ(cmap.validate_intervals(), false);
}

TEST(ChainTest, IntervalMapSanityCheckPartialOverlap) {
    chain::ChainMap cmap;
    // Add intervals with partial overlap
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 300, 50, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 250, 400, 50, true));
    
    EXPECT_EQ(cmap.validate_intervals(), false);
}

TEST(ChainTest, IntervalMapSanityCheckDestOverlap) {
    chain::ChainMap cmap;
    // The dest intervals overlap, but the source intervals do not
    // Dest intervals: chr1_dest: [100, 300), chr1_dest: [100, 200)
    // Source intervals: chr1: [100, 300), chr1: [300, 400)
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 100, 300, 0, true));
    cmap.add_interval("chr1", chain::Interval("chr1_dest", 300, 400, -200, true));
    
    EXPECT_EQ(cmap.validate_intervals(), true);
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
