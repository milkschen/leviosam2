#ifndef LEVIOSAM_MATE_METADATA_HPP
#define LEVIOSAM_MATE_METADATA_HPP

#include <htslib/sam.h>

#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>

namespace LevioSamUtils {

// Only final primary alignment outcomes are indexed. Auxiliary records refer
// to the opposite segment's primary outcome, never another secondary hit.
struct MateOutcome {
    int32_t tid = -1;
    hts_pos_t pos = -1;
    hts_pos_t end = -1;
    uint16_t flag = BAM_FUNMAP;
    uint16_t source_flag = 0;
    uint8_t mapq = 0;
    uint8_t source_orientation = 0;
    bool source_same_reference = false;
};

MateOutcome mate_outcome(const bam1_t* lifted, const bam1_t* source);
void repair_mate_metadata(bam1_t* lifted, const bam1_t* source,
                          const MateOutcome* mate,
                          const MateOutcome* primary = nullptr);

// Temporary BAMs retain the original and final lifted records until every
// primary mate is known. The disk hash index grows on disk, rather than keeping
// one entry per read in RAM. The owner serializes append()/find() with its write
// mutex; replay is single-threaded after all lift workers finish.
class MateMetadataSpool {
   public:
    explicit MateMetadataSpool(uint64_t initial_index_capacity = 65536);
    ~MateMetadataSpool();
    MateMetadataSpool(const MateMetadataSpool&) = delete;
    MateMetadataSpool& operator=(const MateMetadataSpool&) = delete;

    void append(const bam1_t* source, const bam1_t* lifted,
                sam_hdr_t* source_header, sam_hdr_t* lifted_header);
    bool find_mate(const bam1_t* source, MateOutcome& result);
    bool find_primary(const bam1_t* source, MateOutcome& result);
    bool empty() const;
    void begin_replay(sam_hdr_t* source_header, sam_hdr_t* lifted_header);
    bool next(bam1_t* source, bam1_t* lifted,
              sam_hdr_t* source_header, sam_hdr_t* lifted_header);

   private:
    class Impl;
    std::unique_ptr<Impl> impl_;
    uint64_t initial_index_capacity_;
};

}  // namespace LevioSamUtils

#endif
