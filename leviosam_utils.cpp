// #include <thread>
#include "leviosam_utils.hpp"

namespace LevioSamUtils {

void WriteToFastq::init(
    const std::string outpre, const std::string sm, const int mc
) {
    out_fqS.open(outpre + "-deferred-S.fq");
    out_fq1.open(outpre + "-deferred-R1.fq");
    out_fq2.open(outpre + "-deferred-R2.fq");
    split_mode = sm;
    mapq_cutoff = mc;
    r1_db.clear();
    r2_db.clear();
}


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
    // if (verbose > VERBOSE_DEBUG) {
    //     std::cerr << "* new: ";
    //     debug_print_cigar(bam_get_cigar(aln), aln->core.n_cigar);
    // }
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


/* Write a bam1_t object to a FASTQ record. */
void write_fq_from_bam_core(bam1_t* aln, std::ofstream& out_fq){
    out_fq << "@" << bam_get_qname(aln) << "\n";
    out_fq << get_read(aln) << "\n+\n";
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
    out_fq << qual_seq << "\n";
}


/* BAM2FASTQ, with a heuristic algorithm to handle paired-end reads
 * 
 * We don't require the assumption that the dataset is sorted
 * by name, and thus read pairs might not occur in an interleaved
 * format.
 * 
 * We use two unordered_maps to store unmatched R1 and R2 reads.
 * When a segment comes in, we check if the mate is in the map.
 * If so, write both segments to file and mark the mate as 
 * "written"; if not, add the segment to the associated map.
 */
void write_fq_from_bam(bam1_t* aln, WriteToFastq* wfq) {
    if ((aln->core.flag & 256) || // Secondary alignment - no SEQ field
        (aln->core.flag & 512) || // not passing filters
        (aln->core.flag & 1024) || // PCR or optinal duplicate
        (aln->core.flag & 2048)) { // supplementary alignment
        return;
    } 
    std::string n (bam_get_qname(aln));
    if (!(aln->core.flag & 1)) { // e.g. single-end
        write_fq_from_bam_core(aln, wfq->out_fqS);
    } else if (aln->core.flag & 64) { // First segment
        auto find_r2 = wfq->r2_db.find(n);
        if (find_r2 != wfq->r2_db.end()) {
            auto w = find_r2->second.write(wfq->out_fq2, n);
            if (w > 0)
                write_fq_from_bam_core(aln, wfq->out_fq1);
        } else {
            wfq->r1_db.emplace(std::make_pair(n, FastqRecord(aln)));
            // std::cerr << "first+1: " << wfq->r1_db.size() << " " << n << "\n";
        }
    } else if (aln->core.flag & 128) { // Second segment
        auto find_r1 = wfq->r1_db.find(n);
        if (find_r1 != wfq->r1_db.end()) {
            auto w = find_r1->second.write(wfq->out_fq1, n);
            if (w > 0)
                write_fq_from_bam_core(aln, wfq->out_fq2);
        } else {
            wfq->r2_db.emplace(std::make_pair(n, FastqRecord(aln)));
            // std::cerr << "second+1: " << wfq->r2_db.size() << " " << n << "\n";
        }
    } else {
        std::cerr << "Error: issues with read " << n << "\n";
        exit(1);
    }

    if (wfq->r1_db.size() >= 10000) {
        std::cerr << "Clearing R1db (" << wfq->r1_db.size() << ") and R2db (" << wfq->r2_db.size() << ")\n";
        for (auto i: wfq->r1_db) {
            auto n1 = i.first;
            n1 += "/1";
            i.second.write(wfq->out_fqS, n1);
        }
        for (auto i: wfq->r2_db) {
            auto n2 = i.first;
            n2 += "/2";
            i.second.write(wfq->out_fqS, n2);
        }
        wfq->r1_db.clear();
        wfq->r2_db.clear();
    }
}

/* Construct a FastqRecord object */
FastqRecord::FastqRecord(bam1_t* aln) {
    flag = aln->core.flag;
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

/* Write a FastqRecord object to a FASTQ record */
int FastqRecord::write(std::ofstream& out_fq, std::string name) {
    if (written)
        return -1;
    out_fq << "@" << name << "\n";
    out_fq << seq_str << "\n+\n";
    out_fq << qual_str << "\n";
    written = true;
    return 1;
}

};