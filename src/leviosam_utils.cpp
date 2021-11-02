#include <algorithm>
#include <regex>
#include "leviosam_utils.hpp"

namespace LevioSamUtils {

typedef robin_hood::unordered_map<std::string, FastqRecord> fastq_map;

void WriteDeferred::init(
    const std::string outpre, const std::string sm,
    const int mapq, const int isize,
    const float clipped_frac, const int aln_score,
    const std::string of, bam_hdr_t* ihdr
) {
    split_mode = sm;
    min_mapq = mapq;
    max_isize = isize;
    max_clipped_frac = clipped_frac;
    min_aln_score = aln_score;
    hdr_orig = ihdr;
    write_deferred = true;

    std::string out_mode = (of == "sam")? "w" : "wb";
    std::string out_fn = outpre + "-deferred." + of;
    out_fp = sam_open(out_fn.data(), out_mode.data());
    std::string out_fn_orig = outpre + "-unliftable." + of;
    out_fp_orig = sam_open(out_fn_orig.data(), out_mode.data());
    if (sam_hdr_write(out_fp_orig, hdr_orig) < 0) {
        std::cerr << "Error during writing sam_hdr for " << out_fn_orig << "\n";
        exit(1);
    }
}


WriteDeferred::~WriteDeferred() {
    if (write_deferred) {
        sam_close(out_fp);
        sam_close(out_fp_orig);
    }
};


void update_cigar(
    bam1_t* aln, std::vector<uint32_t> &new_cigar
) {
    uint32_t* cigar = bam_get_cigar(aln);
    // Adapted from samtools/sam.c
    // https://github.com/samtools/htslib/blob/2264113e5df1946210828e45d29c605915bd3733/sam.c#L515
    if (aln->core.n_cigar != new_cigar.size()){
        auto cigar_st = (uint8_t*)bam_get_cigar(aln) - aln->data;
        auto fake_bytes = aln->core.n_cigar * 4;
        aln->core.n_cigar = (uint32_t)new_cigar.size();
        auto n_cigar4 = aln->core.n_cigar * 4;
        auto orig_len = aln->l_data;
        if (n_cigar4 > fake_bytes){
            // Check if we need to update `aln->m_data`.
            //
            // Adapted from htslib/sam_internal.h
            // https://github.com/samtools/htslib/blob/31f0a76d338c9bf3a6893b71dd437ef5fcaaea0e/sam_internal.h#L48
            auto new_m_data = (size_t) aln->l_data + n_cigar4 - fake_bytes;
            kroundup32(new_m_data);
            if (new_m_data > aln->m_data){
                auto new_data = static_cast<uint8_t *>(realloc(aln->data, new_m_data));
                if (!new_data){
                    std::cerr << "[Error] Failed to expand a bam1_t struct for " << bam_get_qname(aln) << "\n";
                    std::cerr << "This is likely due to out of memory\n";
                    exit(1);
                }
                aln->data = new_data;
                aln->m_data = new_m_data;
                cigar = bam_get_cigar(aln);
            }
        }
        // Update memory usage of data.
        aln->l_data = aln->l_data - fake_bytes + n_cigar4;
        // Move data to the correct place.
        memmove(aln->data + cigar_st + n_cigar4,
                aln->data + cigar_st + fake_bytes,
                orig_len - (cigar_st + fake_bytes));
        // If new n_cigar is greater, copy the real CIGAR to the right place.
        // Skipped this if new n_cigar is smaller than the original value.
        if (n_cigar4 > fake_bytes){
            memcpy(aln->data + cigar_st,
                   aln->data + (n_cigar4 - fake_bytes) + 8,
                   n_cigar4);
        }
        aln->core.n_cigar = new_cigar.size();
    }
    for (int i = 0; i < aln->core.n_cigar; i++){
        *(cigar + i) = new_cigar[i];
    }
}


void debug_print_cigar(
    uint32_t* cigar, size_t n_cigar
) {
    for (int i = 0; i < n_cigar; i++){
        auto cigar_op_len = bam_cigar_oplen(cigar[i]);
        auto cigar_op = bam_cigar_op(cigar[i]);
        std::cerr << cigar_op_len << bam_cigar_opchr(cigar_op);
    }
    std::cerr << "\n";
}


/* Remove MN:i and MD:z tags from an alignment (bam1_t) object */
void remove_mn_md_tag(bam1_t* aln) {
    uint8_t* ptr = NULL;
    if ((ptr = bam_aux_get(aln, "MD")) != NULL) {
        bam_aux_del(aln, ptr);
    }
    if ((ptr = bam_aux_get(aln, "NM")) != NULL) {
        bam_aux_del(aln, bam_aux_get(aln, "NM"));
    }
}


/* Return the read, reverse complemented if necessary
   Adapted from: https://github.com/samtools/samtools/blob/develop/bam_fastq.c 
*/
static std::string get_read(const bam1_t *rec){
    int len = rec->core.l_qseq + 1;
    char *seq = (char *)bam_get_seq(rec);
    std::string read = "";

    for (int n = 0; n < rec->core.l_qseq; n++) {
        if (rec->core.flag & BAM_FREVERSE)
            read.append(1, seq_nt16_str[seq_comp_table[bam_seqi(seq, n)]]);
        else
            read.append(1, seq_nt16_str[bam_seqi(seq, n)]);
    }
    if (rec->core.flag & BAM_FREVERSE)
        std::reverse(read.begin(), read.end());
    return read;
}


void WriteDeferred::write_deferred_bam(
    bam1_t* aln, bam_hdr_t* hdr
) {
    if ((aln->core.flag & 256) || // Secondary alignment - no SEQ field
        (aln->core.flag & 512) || // not passing filters
        (aln->core.flag & 1024) || // PCR or optinal duplicate
        (aln->core.flag & 2048)) { // supplementary alignment
        return;
    }
    auto w_ret = sam_write1(out_fp, hdr, aln);
    if (w_ret < 0) {
        std::cerr << "[Error] Failed to write record " << bam_get_qname(aln) << "\n";
        exit(1);
    }
}


void WriteDeferred::write_deferred_bam_orig(
    bam1_t* aln
) {
    auto w_ret = sam_write1(out_fp_orig, hdr_orig, aln);
    if (w_ret < 0) {
        std::cerr << "[Error] Failed to write record " << bam_get_qname(aln) << "\n";
        exit(1);
    }
}


/* Construct a FastqRecord object from seq and qual */
FastqRecord::FastqRecord(const std::string& seq, const std::string& qual) {
    seq_str = seq;
    qual_str = qual;
}


/* Construct a FastqRecord object from a BAM record */
FastqRecord::FastqRecord(bam1_t* a) {
    aln = bam_init1();
    if (bam_copy1(aln, a) == NULL) {
        std::cerr << "Failed to copy " << bam_get_qname(a) << "\n";
    };
}


FastqRecord::~FastqRecord() {
    if (aln != NULL) {
        bam_destroy1(aln);
        aln = NULL;
    }
}


/* Write a FastqRecord object to a FASTQ file */
int FastqRecord::write(std::ofstream& out_fq, std::string name) {
    if (aln != NULL) {
        seq_str = get_read(aln);
        std::string qual_seq("");
        uint8_t* qual = bam_get_qual(aln);
        if (qual[0] == 255) qual_seq = "*";
        else {
            for (auto i = 0; i < aln->core.l_qseq; ++i) {
                qual_seq += (char) (qual[i] + 33);
            }
        }
        if (aln->core.flag & BAM_FREVERSE)
            std::reverse(qual_seq.begin(), qual_seq.end());
        qual_str = qual_seq;
    }
    out_fq << "@" << name << "\n";
    out_fq << seq_str << "\n+\n";
    out_fq << qual_str << "\n";
    return 1;
}


/* Read a SAM/BAM file, write properly paired alignements to a BAM file
 * and return the rest as a fastq_map
 */
fastq_map read_deferred_bam(
    samFile* dsam_fp, samFile* out_dsam_fp, bam_hdr_t* hdr,
    std::ofstream& out_r1_fp, std::ofstream& out_r2_fp
) {
    fastq_map reads1, reads2;
    bam1_t* aln = bam_init1();

    while (sam_read1(dsam_fp, hdr, aln) > 0) {
        bam1_core_t c = aln->core;
        // The following categories of reads are excluded by this method
        if ((c.flag & 256) || // Secondary alignment - no SEQ field
            (c.flag & 512) || // not passing filters
            (c.flag & 1024) || // PCR or optinal duplicate
            (c.flag & 2048)) { // supplementary alignment
            continue;
        }
        std::string qname = bam_get_qname(aln);
        LevioSamUtils::FastqRecord fq = LevioSamUtils::FastqRecord(aln);
        if (c.flag & 64) { // first segment
            auto search = reads2.find(qname);
            if (search != reads2.end()) {
                if (sam_write1(out_dsam_fp, hdr, search->second.aln) < 0  ||
                    sam_write1(out_dsam_fp, hdr, aln) < 0) {
                    std::cerr << "[Error] Failed to write record " << qname << "\n";
                    exit(1);
                }
                fq.write(out_r1_fp, qname);
                search->second.write(out_r2_fp, qname);
                reads2.erase(search);
            } else {
                reads1.emplace(std::make_pair(qname, fq));
            }
        } else if (c.flag & 128) { // second segment
            auto search = reads1.find(qname);
            if (search != reads1.end()) {
                if (sam_write1(out_dsam_fp, hdr, aln) < 0 ||
                    sam_write1(out_dsam_fp, hdr, search->second.aln) < 0) {
                    std::cerr << "[Error] Failed to write record " << qname << "\n";
                    exit(1);
                }
                search->second.write(out_r1_fp, qname);
                fq.write(out_r2_fp, qname);
                reads1.erase(search);
            } else {
                reads2.emplace(std::make_pair(qname, fq));
            }
        } else {
            std::cerr << "[Error] Alignment " << qname << " is not paired.\n";
            exit(1);
        }
    }
    std::cerr << "Num. unpaired reads1: " << reads1.size() << "\n";
    std::cerr << "Num. unpaired reads2: " << reads2.size() << "\n";
    for (auto& r: reads2) {
        reads1.emplace(std::make_pair(r.first, r.second.aln));
    }
    std::cerr << "Num. merged unpaired: " << reads1.size() << "\n";
    return reads1;
}


/* Read a single-end FASTQ file and return a fastq_map */
fastq_map read_unpaired_fq(const std::string& fq_fname) {
    fastq_map reads;
    std::ifstream fastq_fp(fq_fname);
    std::string line;
    int i = 0;
    std::string name, seq;
    while (getline (fastq_fp, line)) {
        if (i % 4 == 0) {
            name = line.substr(1);
        } else if (i % 4 == 1) {
            seq = line;
        } else if (i % 4 == 3) {
            reads.emplace(
                std::make_pair(name,
                               LevioSamUtils::FastqRecord(seq, line)));
            name = "";
            seq = "";
        }
        i++;
    }
    std::cerr << "Number of singletons: " << reads.size() << "\n";
    fastq_fp.close();
    return reads;
}


/* Split a string on a delimiter
 * From: https://stackoverflow.com/a/64886763
 */
std::vector<std::string> split_str(
    const std::string str, const std::string regex_str
) {
    std::regex regexz(regex_str);
    std::vector<std::string> list(
        std::sregex_token_iterator(str.begin(), str.end(), regexz, -1),
        std::sregex_token_iterator());
    return list;
}


/* Update SEQ and QUAL if in a reversed chain */
int reverse_seq_and_qual(bam1_t* aln) {
    bam1_core_t* c = &(aln->core);
    uint8_t *seq = bam_get_seq(aln);
    uint8_t *qual = bam_get_qual(aln);
    std::vector<int> read_nt16;
    std::vector<uint8_t> qual_tmp;
    for (int n = 0; n < c->l_qseq; n++) {
        read_nt16.push_back(seq_comp_table[bam_seqi(seq, n)]);
        qual_tmp.push_back(qual[n]);
    }
    // Make sure the length of the reversed seq is concordant with the original seq
    if (read_nt16.size() != c->l_qseq) {
        return -1;
    }
    for (int n = 0; n < c->l_qseq; n++) {
        bam_set_seqi(seq, n, read_nt16[c->l_qseq - n - 1]);
        qual[n] = qual_tmp[c->l_qseq - n - 1];
    }

    return 0;
}

};
