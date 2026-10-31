#include "mate_metadata.hpp"

#include <htslib/hfile.h>

#include <algorithm>
#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <unistd.h>

#include "leviosam_utils.hpp"

namespace LevioSamUtils {
namespace {

[[noreturn]] void fail(const std::string& operation) {
    throw std::runtime_error("mate metadata: " + operation + ": " +
                             std::strerror(errno));
}

class TempFile {
   public:
    TempFile() : file_(nullptr) {
        const char* configured = std::getenv("TMPDIR");
        std::string directory = configured && *configured ? configured : "/tmp";
        std::string path = directory + "/leviosam2-mates-XXXXXX";
        // mkstemp creates mode 0600. Unlink immediately: failures, exceptions
        // and process termination cannot leave genome/read data behind.
        int fd = mkstemp(&path[0]);
        if (fd < 0) fail("cannot create private temporary file");
        if (unlink(path.c_str()) != 0) {
            close(fd);
            fail("cannot unlink private temporary file");
        }
        file_ = fdopen(fd, "w+b");
        if (!file_) {
            close(fd);
            fail("cannot open private temporary file");
        }
    }
    ~TempFile() { if (file_) std::fclose(file_); }
    FILE* get() const { return file_; }
    void swap(TempFile& other) { std::swap(file_, other.file_); }
    void seek(uint64_t offset) {
        if (offset > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
            fseeko(file_, static_cast<off_t>(offset), SEEK_SET) != 0)
            fail("temporary-file seek failed");
    }
    uint64_t end() {
        if (fseeko(file_, 0, SEEK_END) != 0) fail("temporary-file seek failed");
        off_t offset = ftello(file_);
        if (offset < 0) fail("temporary-file position failed");
        return static_cast<uint64_t>(offset);
    }
    void read(void* data, size_t size) {
        if (std::fread(data, 1, size, file_) != size)
            fail("temporary-file read failed");
    }
    void write(const void* data, size_t size) {
        if (std::fwrite(data, 1, size, file_) != size)
            fail("temporary-file write failed");
    }
    void resize(uint64_t size) {
        if (size > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) ||
            ftruncate(fileno(file_), static_cast<off_t>(size)) != 0)
            fail("temporary index allocation failed");
    }
   private:
    FILE* file_;
};

samFile* bam_handle(TempFile& file, const char* mode) {
    // A duplicate descriptor gives HTSlib independent ownership while tmpfile
    // keeps the unlinked file alive between the write and replay passes.
    int fd = dup(fileno(file.get()));
    if (fd < 0) fail("cannot duplicate temporary BAM descriptor");
    hFILE* handle = hdopen(fd, mode);
    if (!handle) {
        close(fd);
        fail("cannot open temporary BAM descriptor");
    }
    samFile* sam = hts_hopen(handle, "mate-metadata-spool.bam", mode);
    if (!sam) {
        int close_status = hclose(handle);
        if (close_status < 0) fail("cannot close failed temporary BAM descriptor");
        fail("cannot open temporary BAM");
    }
    return sam;
}

std::string key(const bam1_t* source, bool opposite) {
    const bool first = source->core.flag & BAM_FREAD1;
    const bool second = source->core.flag & BAM_FREAD2;
    if (first == second)
        throw std::runtime_error("mate metadata: paired record " +
            std::string(bam_get_qname(source)) +
            " must identify exactly one of READ1 and READ2");
    std::string result(bam_get_qname(source));
    result.push_back('\0');
    uint8_t* rg = bam_aux_get(source, "RG");
    if (rg) {
        const char* group = bam_aux2Z(rg);
        if (!group) throw std::runtime_error("mate metadata: RG must be a string");
        result += group;
    }
    result.push_back('\0');
    result.push_back((first != opposite) ? '1' : '2');
    return result;
}

uint64_t hash_key(const std::string& text) {
    // FNV-1a, followed by an avalanche so low table bits are well distributed.
    uint64_t hash = UINT64_C(14695981039346656037);
    for (unsigned char ch : text) {
        hash ^= ch;
        hash *= UINT64_C(1099511628211);
    }
    hash ^= hash >> 33;
    hash *= UINT64_C(0xff51afd7ed558ccd);
    hash ^= hash >> 33;
    return hash;
}

uint8_t orientation(bool reverse, bool mate_reverse, bool left) {
    if (reverse == mate_reverse) return 0;  // same-strand library
    return (left ? reverse : mate_reverse) ? 2 : 1;  // RF / FR
}

bool proper_pair(const MateOutcome& primary, const MateOutcome& mate) {
    if ((primary.flag & BAM_FUNMAP) || (mate.flag & BAM_FUNMAP) ||
        primary.tid != mate.tid || !primary.source_same_reference ||
        !(primary.source_flag & BAM_FPROPER_PAIR) ||
        !(mate.source_flag & BAM_FPROPER_PAIR)) return false;
    bool left = primary.pos < mate.pos ||
                (primary.pos == mate.pos && (primary.flag & BAM_FREAD1));
    return primary.source_orientation == orientation(primary.flag & BAM_FREVERSE,
                                                       mate.flag & BAM_FREVERSE, left);
}

struct IndexSlot {
    uint64_t hash = 0;
    uint64_t offset = 0;  // record offset + 1; zero means empty
};

class OutcomeIndex {
   public:
    explicit OutcomeIndex(uint64_t capacity) : capacity_(capacity) {
        index_.resize(capacity_ * sizeof(IndexSlot));
    }
    bool find(const std::string& key, MateOutcome& result) {
        IndexSlot slot = locate(key, hash_key(key)).second;
        if (!slot.offset) return false;
        read_record(slot.offset, nullptr, &result);
        return true;
    }
    void insert(const std::string& key, const MateOutcome& outcome) {
        if ((entries_ + 1) * 10 > capacity_ * 7) grow();
        uint64_t hash = hash_key(key);
        auto found = locate(key, hash);
        if (found.second.offset)
            throw std::runtime_error("mate metadata: ambiguous duplicate primary "
                "alignment for QNAME/RG/read end: " + key.substr(0, key.find('\0')));
        if (key.size() > std::numeric_limits<uint32_t>::max())
            throw std::runtime_error("mate metadata: primary identity is too long");
        uint32_t size = static_cast<uint32_t>(key.size());
        uint64_t offset = records_.end();
        records_.write(&size, sizeof(size));
        records_.write(key.data(), key.size());
        records_.write(&outcome, sizeof(outcome));
        IndexSlot slot{hash, offset + 1};
        index_.seek(found.first * sizeof(IndexSlot));
        index_.write(&slot, sizeof(slot));
        ++entries_;
    }
   private:
    TempFile index_, records_;
    uint64_t capacity_;
    uint64_t entries_ = 0;

    void read_record(uint64_t offset, std::string* key, MateOutcome* outcome) {
        records_.seek(offset - 1);
        uint32_t size;
        records_.read(&size, sizeof(size));
        if (key) {
            key->resize(size);
            records_.read(&(*key)[0], size);
        } else {
            records_.seek(offset - 1 + sizeof(size) + size);
        }
        if (outcome) records_.read(outcome, sizeof(*outcome));
    }
    std::pair<uint64_t, IndexSlot> locate(const std::string& key, uint64_t hash) {
        uint64_t i = hash & (capacity_ - 1);
        while (true) {
            IndexSlot slot;
            index_.seek(i * sizeof(slot));
            index_.read(&slot, sizeof(slot));
            if (!slot.offset) return {i, slot};
            if (slot.hash == hash) {
                std::string candidate;
                read_record(slot.offset, &candidate, nullptr);
                if (candidate == key) return {i, slot};
            }
            i = (i + 1) & (capacity_ - 1);
        }
    }
    void grow() {
        if (capacity_ > std::numeric_limits<uint64_t>::max() / (2 * sizeof(IndexSlot)))
            throw std::runtime_error("mate metadata: temporary index is too large");
        uint64_t next_capacity = capacity_ * 2;
        TempFile next;
        next.resize(next_capacity * sizeof(IndexSlot));
        for (uint64_t i = 0; i < capacity_; ++i) {
            IndexSlot slot;
            index_.seek(i * sizeof(slot));
            index_.read(&slot, sizeof(slot));
            if (!slot.offset) continue;
            uint64_t j = slot.hash & (next_capacity - 1);
            while (true) {
                IndexSlot candidate;
                next.seek(j * sizeof(candidate));
                next.read(&candidate, sizeof(candidate));
                if (!candidate.offset) break;
                j = (j + 1) & (next_capacity - 1);
            }
            next.seek(j * sizeof(slot));
            next.write(&slot, sizeof(slot));
        }
        index_.swap(next);
        capacity_ = next_capacity;
    }
};

}  // namespace

MateOutcome mate_outcome(const bam1_t* lifted, const bam1_t* source) {
    MateOutcome result;
    result.tid = lifted->core.tid;
    result.pos = lifted->core.pos;
    result.end = bam_endpos(lifted);
    result.flag = lifted->core.flag;
    result.source_flag = source->core.flag;
    result.mapq = lifted->core.qual;
    bool left = source->core.pos < source->core.mpos ||
        (source->core.pos == source->core.mpos && (source->core.flag & BAM_FREAD1));
    result.source_orientation = orientation(source->core.flag & BAM_FREVERSE,
                                             source->core.flag & BAM_FMREVERSE, left);
    result.source_same_reference = source->core.tid == source->core.mtid;
    return result;
}

void repair_mate_metadata(bam1_t* lifted, const bam1_t* source,
                          const MateOutcome* mate, const MateOutcome* primary) {
    bam1_core_t& c = lifted->core;
    if (!(c.flag & BAM_FPAIRED)) return;
    bool mapped = !(c.flag & BAM_FUNMAP);
    // Proper pairing is reassessed from both final alignments below, not an
    // earlier independent estimate that may have marked the mate unliftable.
    c.flag &= ~BAM_FPROPER_PAIR;
    c.isize = 0;
    if (!mate) {
        // Partial BAM with neither a primary mate nor an exact MC fallback.
        // The source knows mappedness but not the mate's target location.
        c.flag &= ~(BAM_FMUNMAP | BAM_FMREVERSE);
        c.flag |= source->core.flag & BAM_FMUNMAP;
        c.mtid = -1;
        c.mpos = -1;
        if (c.flag & BAM_FMUNMAP) {
            c.mtid = mapped ? c.tid : -1;
            c.mpos = mapped ? c.pos : -1;
        }
        if (!mapped) c.tid = c.pos = -1;
        remove_aux_tag(lifted, "MQ");
        const std::string status = std::string(mapped ? "L" : "UM") + "_" +
                                   ((c.flag & BAM_FMUNMAP) ? "UM" : "L");
        remove_aux_tag(lifted, "LO");
        if (bam_aux_update_str(lifted, "LO", status.size() + 1, status.c_str()) < 0)
            throw std::runtime_error("mate metadata: cannot update LO tag");
        return;
    }
    bool mate_mapped = !(mate->flag & BAM_FUNMAP);
    c.flag &= ~(BAM_FMUNMAP | BAM_FMREVERSE);
    if (!mate_mapped) c.flag |= BAM_FMUNMAP;
    if (mate_mapped && (mate->flag & BAM_FREVERSE)) c.flag |= BAM_FMREVERSE;
    c.mtid = mate_mapped ? mate->tid : (mapped ? c.tid : -1);
    c.mpos = mate_mapped ? mate->pos : (mapped ? c.pos : -1);
    if (!mapped) {
        c.tid = mate_mapped ? mate->tid : -1;
        c.pos = mate_mapped ? mate->pos : -1;
    }
    bool nonprimary = c.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY);
    // Supplementary/secondary records may inherit the primary pair's proper
    // status even at a different locus (e.g. BWA). Judge that status using the
    // two actual primaries, not the supplementary record's own chromosome.
    MateOutcome own = mate_outcome(lifted, source);
    const MateOutcome* pair_primary = nonprimary ? primary : &own;
    if (mapped && pair_primary &&
        (!nonprimary || (source->core.flag & BAM_FPROPER_PAIR)) &&
        proper_pair(*pair_primary, *mate))
        c.flag |= BAM_FPROPER_PAIR;
    if (mapped && mate_mapped && c.tid == mate->tid) {
        bool left = c.pos < mate->pos ||
                    (c.pos == mate->pos && (c.flag & BAM_FREAD1));
        if (!nonprimary) {
            hts_pos_t span = std::max(bam_endpos(lifted), mate->end) -
                             std::min(c.pos, mate->pos);
            c.isize = left ? span : -span;
        }
    }
    if (bam_aux_get(lifted, "MQ")) {
        if (bam_aux_update_int(lifted, "MQ", mate_mapped ? mate->mapq : 0) < 0)
            throw std::runtime_error("mate metadata: cannot update MQ tag");
    }
    const std::string status = std::string(mapped ? "L" : "UM") + "_" +
                               (mate_mapped ? "L" : "UM");
    remove_aux_tag(lifted, "LO");
    if (bam_aux_update_str(lifted, "LO", status.size() + 1, status.c_str()) < 0)
        throw std::runtime_error("mate metadata: cannot update LO tag");
}

class MateMetadataSpool::Impl {
   public:
    explicit Impl(uint64_t capacity) : outcomes(capacity) {}
    TempFile source_file, lifted_file;
    OutcomeIndex outcomes;
    samFile* source = nullptr;
    samFile* lifted = nullptr;
    uint64_t count = 0;
    bool replay = false;
    uint64_t replayed = 0;
    ~Impl() {
        if (source) sam_close(source);
        if (lifted) sam_close(lifted);
    }
};

MateMetadataSpool::MateMetadataSpool(uint64_t capacity)
    : initial_index_capacity_(capacity) {
    if (capacity < 4 || (capacity & (capacity - 1)) != 0 ||
        capacity > static_cast<uint64_t>(std::numeric_limits<off_t>::max()) / sizeof(IndexSlot))
        throw std::runtime_error("mate metadata: invalid temporary index capacity");
}
MateMetadataSpool::~MateMetadataSpool() = default;
bool MateMetadataSpool::empty() const { return !impl_ || !impl_->count; }

void MateMetadataSpool::append(const bam1_t* source, const bam1_t* lifted,
                               sam_hdr_t* source_header, sam_hdr_t* lifted_header) {
    if (!impl_) {
        impl_.reset(new Impl(initial_index_capacity_));
        std::cerr << "[I::lift_run] Reconciling paired alignments through private "
                     "temporary BAMs and a disk-backed primary-mate index\n";
        impl_->source = bam_handle(impl_->source_file, "wb1");
        impl_->lifted = bam_handle(impl_->lifted_file, "wb1");
        if (sam_hdr_write(impl_->source, source_header) < 0 ||
            sam_hdr_write(impl_->lifted, lifted_header) < 0)
            fail("cannot write temporary BAM header");
    }
    if (!(source->core.flag & (BAM_FSECONDARY | BAM_FSUPPLEMENTARY)))
        impl_->outcomes.insert(key(source, false), mate_outcome(lifted, source));
    if (sam_write1(impl_->source, source_header, source) < 0 ||
        sam_write1(impl_->lifted, lifted_header, lifted) < 0)
        fail("cannot write temporary BAM record");
    ++impl_->count;
}

bool MateMetadataSpool::find_mate(const bam1_t* source, MateOutcome& result) {
    return impl_ && impl_->outcomes.find(key(source, true), result);
}

bool MateMetadataSpool::find_primary(const bam1_t* source, MateOutcome& result) {
    return impl_ && impl_->outcomes.find(key(source, false), result);
}

void MateMetadataSpool::begin_replay(sam_hdr_t*, sam_hdr_t*) {
    if (empty()) return;
    samFile* source = impl_->source;
    samFile* lifted = impl_->lifted;
    impl_->source = impl_->lifted = nullptr;
    int source_status = sam_close(source), lifted_status = sam_close(lifted);
    if (source_status < 0 || lifted_status < 0)
        fail("cannot finish temporary BAM");
    impl_->source_file.seek(0);
    impl_->lifted_file.seek(0);
    impl_->source = bam_handle(impl_->source_file, "r");
    impl_->lifted = bam_handle(impl_->lifted_file, "r");
    sam_hdr_t* sh = sam_hdr_read(impl_->source);
    sam_hdr_t* lh = sam_hdr_read(impl_->lifted);
    bool valid = sh && lh;
    sam_hdr_destroy(sh);
    sam_hdr_destroy(lh);
    if (!valid) fail("cannot read temporary BAM header");
    impl_->replay = true;
}

bool MateMetadataSpool::next(bam1_t* source, bam1_t* lifted,
                             sam_hdr_t* source_header, sam_hdr_t* lifted_header) {
    if (empty() || !impl_->replay) return false;
    int sr = sam_read1(impl_->source, source_header, source);
    int lr = sam_read1(impl_->lifted, lifted_header, lifted);
    if (sr == -1 && lr == -1) {
        if (impl_->replayed != impl_->count)
            throw std::runtime_error("mate metadata: truncated temporary BAMs");
        return false;
    }
    if (sr < 0 || lr < 0) fail("cannot read matching temporary BAM records");
    ++impl_->replayed;
    return true;
}

}  // namespace LevioSamUtils
