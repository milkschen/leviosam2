/*
 * liftover.hpp
 *
 * Classes and routines for translating (lifting over) coordinates between
 * two aligned sequences
 *
 * Authors: Taher Mun, Nae-Chyun Chen, Ben Langmead
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#ifndef LIFTOVER_HPP
#define LIFTOVER_HPP

#include <htslib/sam.h>
#include <htslib/vcf.h>

#include <cstdio>
#include <iostream>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
#include <thread>

#include "aln.hpp"
#include "bed.hpp"
#include "chain.hpp"
#include "leviosam_utils.hpp"
#include "version.hpp"
#include "yaml.hpp"

#define VERBOSE_CLEAN 0
#define VERBOSE_INFO 1
#define VERBOSE_DEBUG 2
#define VERBOSE_DEV 3

using NameMap = std::vector<std::pair<std::string, std::string>>;
using LengthMap = std::vector<std::pair<std::string, int32_t>>;
using BedMap =
    robin_hood::unordered_map<std::string, IITree<std::size_t, bool>>;

static inline void die(std::string msg) {
    fprintf(stderr, "%s\n", msg.data());
    exit(1);
}

extern "C" {
int bam_fillmd1(bam1_t *b, const char *ref, int flag, int quiet_mode);
}

struct lift_opts {
    std::string vcf_fname = "";
    std::string chain_fname = "";
    std::string chainmap_fname = "";
    std::string sample = "";
    std::string outpre = "";
    std::string out_format = "sam";
    std::string lift_fname = "";
    std::string sam_fname = "";
    std::string cmd = "";
    std::string haplotype = "0";
    std::string split_mode = "";
    std::vector<std::pair<std::string, float>> split_rules;
    int allowed_cigar_changes = 0;
    int threads = 1;
    int chunk_size = 256;
    int verbose = 0;
    int md_flag = 0;
    std::string ref_name = "";
    std::string realign_yaml = "";
    Aln::AlnOpts aln_opts;
    NameMap name_map;
    LengthMap length_map;
    BedUtils::Bed bed_defer_source;
    BedUtils::Bed bed_defer_dest;
    BedUtils::Bed bed_commit_source;
    BedUtils::Bed bed_commit_dest;
    // 0: take any overlap, >1: base pairs, 1>=value>0: fraction
    float bed_isec_threshold = 0;  
};

#define LIFT_R_L 0           // unpaired, liftable
#define LIFT_R_UL 1          // unpaired, un-liftable
#define LIFT_R_UM 2          // unpaired, unmapped
#define LIFT_R1_L_R2_L 3     // paired, R1 liftable, R2 liftable
#define LIFT_R1_L_R2_UL 4    // paired, R1 liftable, R2 un-liftable
#define LIFT_R1_L_R2_UM 5    // paired, R1 liftable, R2 unmapped
#define LIFT_R1_UL_R2_L 6    // paired, R1 un-liftable, R2 liftable
#define LIFT_R1_UL_R2_UL 7   // paired, R1 un-liftable, R2 un-liftable
#define LIFT_R1_UL_R2_UM 8   // paired, R1 un-liftable, R2 unmapped
#define LIFT_R1_UM_R2_L 9    // paired, R1 unmapped, R2 liftable
#define LIFT_R1_UM_R2_UL 10  // paired, R1 unmapped, R2 un-liftable
#define LIFT_R1_UM_R2_UM 11  // paired, R1 unmapped, R2 unmapped

void print_main_help_msg();
void print_lift_help_msg();
void print_serialize_help_msg();

bool check_split_rule(std::string rule);
bool add_split_rule(std::vector<std::pair<std::string, float>> &split_rules,
                    std::string s);

namespace lift {
// Serialization
class Name2IdxMap : public std::unordered_map<std::string, int> {
   public:
    size_t serialize(std::ofstream &out) {
        size_t size = 0;
        size_t map_size = this->size();
        out.write(reinterpret_cast<char *>(&map_size), sizeof(map_size));
        for (auto s : *this) {
            size_t str_size = s.first.size();
            out.write(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            out.write(reinterpret_cast<const char *>(s.first.data()), str_size);
            out.write(reinterpret_cast<char *>(&s.second), sizeof(s.second));
            size += sizeof(str_size) + str_size + sizeof(s.second);
        }
        return size;
    }

    void load(std::ifstream &in) {
        size_t map_size;
        in.read(reinterpret_cast<char *>(&map_size), sizeof(map_size));
        for (auto i = 0; i < map_size; ++i) {
            size_t str_size;
            in.read(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            std::vector<char> buf(str_size);
            in.read(reinterpret_cast<char *>(buf.data()), str_size);
            key_type key(buf.begin(), buf.end());
            mapped_type value;
            in.read(reinterpret_cast<char *>(&value), sizeof(value));
            // fprintf(stderr, "%s\t%d\n", key.data(), value);
            (*this)[key] = value;
        }
    }
};

class Name2NameMap : public std::unordered_map<std::string, std::string> {
   public:
    size_t serialize(std::ofstream &out) {
        size_t size = 0;
        size_t map_size = this->size();
        out.write(reinterpret_cast<char *>(&map_size), sizeof(map_size));
        for (auto s : *this) {
            size_t str_size = s.first.size();
            out.write(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            out.write(reinterpret_cast<const char *>(s.first.data()), str_size);
            size += sizeof(str_size) + str_size;
            str_size = s.second.size();
            out.write(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            out.write(reinterpret_cast<const char *>(s.second.data()),
                      str_size);
            size += sizeof(str_size) + str_size;
        }
        return size;
    }

    void load(std::ifstream &in) {
        size_t map_size;
        in.read(reinterpret_cast<char *>(&map_size), sizeof(map_size));
        for (auto i = 0; i < map_size; ++i) {
            size_t str_size;
            in.read(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            std::vector<char> buf(str_size);
            in.read(reinterpret_cast<char *>(buf.data()), str_size);
            key_type key(buf.begin(), buf.end());
            in.read(reinterpret_cast<char *>(&str_size), sizeof(str_size));
            buf.clear();
            buf.resize(str_size);
            in.read(reinterpret_cast<char *>(buf.data()), str_size);
            mapped_type value(buf.begin(), buf.end());
            (*this)[key] = value;
        }
    }
};

// liftover data structure for a single sequence
class Lift {
   public:
    Lift() { init_rs_sls(); }

    Lift(std::ifstream &in) { this->load(in); }

    // construct liftover from bitvectors
    // i: contains 1 at each insertion in s2 wrt s1
    // d: contains 1 at each deletion  ""
    // s: contains 1 at each mismatch in alignment bwtn s1 and s1
    Lift(sdsl::bit_vector i, sdsl::bit_vector d, sdsl::bit_vector s)
        : ins(i), del(d), snp(s) {
        // create rank and select data-structures for everything
        init_rs_sls();
    }

    Lift(const Lift &rhs) : ins(rhs.ins), del(rhs.del), snp(rhs.snp) {
        init_rs_sls();
    }

    Lift(Lift &&rhs)
        : ins(std::move(rhs.ins)),
          del(std::move(rhs.del)),
          snp(std::move(rhs.snp)) {
        init_rs_sls();
    }

    Lift &operator=(const Lift &rhs) {
        ins = rhs.ins;
        del = rhs.del;
        snp = rhs.snp;
        init_rs_sls();
        return *this;
    }

    Lift &operator=(Lift &&rhs) {
        ins = std::move(rhs.ins);
        del = std::move(rhs.del);
        snp = std::move(rhs.snp);
        init_rs_sls();
        return *this;
    }

    // translate coordinate in s2-space to s1-space
    size_t lift_pos(size_t p) const {  // s2_to_s1
        return ins_rs0(del_sls0(p + 1));
    }

    // translate coordinate in s1-space to s2-space
    // size_t s1_to_s2(size_t p) const {
    //     return del_rs0(ins_sls0(p+1));
    // }

    /* Core function to lift a CIGAR string from s1 to s2 space.
     *
     * Input:
     *   An alignment struct we'd like to update.
     *   This function is not going to update the bam1_t stuct.
     *
     * Output:
     *   A vector, each element is an unit32_t corresponding to htslib CIGAR
     * value
     */
    std::vector<uint32_t> lift_cigar_core(bam1_t *b) {
        auto x = del_sls0(b->core.pos + 1);
        int y = 0;  // read2alt cigar pointer
        uint32_t *cigar = bam_get_cigar(b);
        std::vector<uint32_t> cigar_ops;
        std::vector<uint32_t> new_cigar_ops;
        for (int i = 0; i < b->core.n_cigar; ++i) {
            for (int j = 0; j < bam_cigar_oplen(cigar[i]); ++j) {
                cigar_ops.push_back(bam_cigar_op(cigar[i]));
            }
        }
        int iters = cigar_ops.size();
        while (y < iters) {
            int cop = cigar_ops[y];
            if (del[x]) {  // skip ahead in alt2ref cigar
                if (new_cigar_ops.empty()) {
                    new_cigar_ops.push_back(
                        BAM_CDEL);  // TODO: maybe change this to a padded skip?
                } else {
                    if (new_cigar_ops.back() == BAM_CINS) {
                        new_cigar_ops[new_cigar_ops.size() - 1] = BAM_CMATCH;
                    } else if (new_cigar_ops.back() == BAM_CREF_SKIP) {
                        new_cigar_ops.push_back(BAM_CREF_SKIP);
                    } else {
                        new_cigar_ops.push_back(BAM_CDEL);
                    }
                }
                ++x;
            } else if (cop == BAM_CINS) {  // skip ahead in read2alt cigar
                if (!new_cigar_ops.empty()) {
                    if (new_cigar_ops.back() == BAM_CDEL) {
                        new_cigar_ops[new_cigar_ops.size() - 1] = BAM_CMATCH;
                    } else {
                        new_cigar_ops.push_back(BAM_CINS);
                    }
                } else {
                    new_cigar_ops.push_back(BAM_CINS);
                }
                ++y;
            } else if (cop == BAM_CSOFT_CLIP || cop == BAM_CHARD_CLIP ||
                       cop == BAM_CPAD) {
                // TODO: check the case where S is followed by D
                new_cigar_ops.push_back(cop);
                ++y;
            } else if (cop == BAM_CBACK) {
                // skip. We don't support B Ops.
                fprintf(stderr, "Warning: B operators are not supported\n");
                ++y;
            } else if (cop == BAM_CMATCH || cop == BAM_CDIFF ||
                       cop == BAM_CEQUAL) {  // M
                if (ins[x]) {
                    if (new_cigar_ops.empty()) {
                        new_cigar_ops.push_back(BAM_CINS);
                    } else if (new_cigar_ops.back() == BAM_CDEL) {
                        new_cigar_ops[new_cigar_ops.size() - 1] = BAM_CMATCH;
                    } else {
                        new_cigar_ops.push_back(BAM_CINS);
                    }
                } else {
                    new_cigar_ops.push_back(BAM_CMATCH);
                }
                ++x;
                ++y;
            } else if (cop == BAM_CREF_SKIP) {  // N
                // we separate this from the D branch in case we have to do
                // something special ie. spliced alignments
                if (!ins[x]) {
                    if (new_cigar_ops.empty()) {
                        new_cigar_ops.push_back(BAM_CREF_SKIP);
                    } else if (new_cigar_ops.back() == BAM_CINS) {
                        new_cigar_ops[new_cigar_ops.size() - 1] = BAM_CMATCH;
                    } else {
                        new_cigar_ops.push_back(BAM_CREF_SKIP);
                    }
                }
                ++x;
                ++y;
            } else if (cop == BAM_CDEL) {  // D
                if (!ins[x]) {
                    if (new_cigar_ops.empty()) {
                        new_cigar_ops.push_back(BAM_CDEL);
                    } else if (new_cigar_ops.back() == BAM_CINS) {
                        new_cigar_ops[new_cigar_ops.size() - 1] = BAM_CMATCH;
                    } else {
                        new_cigar_ops.push_back(BAM_CDEL);
                    }
                }
                ++x;
                ++y;
            }
        }
        return new_cigar_ops;
    }

    /* Convert a CIGAR string in an s2 alignment to the corresponding CIGAR
     * string in the s1 alignment
     *
     * Input: bam1_t alignment record against s2 sequence
     *
     * Output:
     *   The bam1_t struct is updated if there's a change in its CIGAR.
     *   bam1_t->core.n_cigar specifies numbers of operators in a CIGAR;
     *   bam_get_cigar(bam1_t) points to the occurrences.
     */
    void lift_cigar(bam1_t *b) {
        uint32_t *cigar = bam_get_cigar(b);
        // Get lifted CIGAR.
        auto new_cigar_ops = lift_cigar_core(b);
        // Prepare to update CIGAR in the bam1_t struct.
        std::vector<int> cigar_op_len;
        std::vector<char> cigar_op;
        int z = 1;
        for (int i = 1; i < new_cigar_ops.size(); i++) {
            if (new_cigar_ops[i] == new_cigar_ops[i - 1]) {
                z++;
            } else {
                cigar_op_len.push_back(z);
                cigar_op.push_back(new_cigar_ops[i - 1]);
                z = 1;
            }
        }
        cigar_op_len.push_back(z);
        cigar_op.push_back(new_cigar_ops[new_cigar_ops.size() - 1]);
        std::vector<uint32_t> new_cigar;
        for (int i = 0; i < cigar_op_len.size(); i++) {
            new_cigar.push_back(bam_cigar_gen(cigar_op_len[i], cigar_op[i]));
        }

        LevioSamUtils::update_cigar(b, new_cigar);
    }

    // returns size of s1 sequence
    size_t s1_len() {
        return ins[ins.size() - 1] ? ins_rs0(ins.size() - 1)
                                   : ins_rs0(ins.size() - 1) + 1;
    }

    // returns size of s2 sequence
    // size_t s2_len() {
    //     return del[del.size() - 1] ? del_rs0(del.size() - 1) :
    //     del_rs0(del.size() - 1) + 1;
    // }

    // Save to stream
    size_t serialize(std::ofstream &out) const {
        size_t size = 0;
        size += ins.serialize(out);
        size += del.serialize(out);
        size += snp.serialize(out);
        return size;
    }

    // Load from stream
    void load(std::istream &in) {
        ins.load(in);
        del.load(in);
        snp.load(in);
        init_rs_sls();
    }

   private:
    void init_rs_sls() {
        sdsl::util::init_support(ins_rs0, &ins);
        sdsl::util::init_support(snp_rs0, &snp);
        sdsl::util::init_support(del_sls0, &del);
        sdsl::util::init_support(snp_sls0, &snp);
    }

    sdsl::sd_vector<> ins;
    sdsl::sd_vector<> del;
    sdsl::sd_vector<> snp;
    sdsl::sd_vector<>::rank_0_type ins_rs0;
    sdsl::sd_vector<>::rank_0_type snp_rs0;
    sdsl::sd_vector<>::select_0_type del_sls0;
    sdsl::sd_vector<>::select_0_type snp_sls0;
};

// for containing and handling 1+ Lift objects
// use this to access Lifts by sequence name
// currently supports construction from a vcf file, other construction methods
// to be supported later.
class LiftMap {
   public:
    LiftMap() {}

    ~LiftMap() {}

    LiftMap(std::ifstream &in) { this->load(in); }

    // copy constructor
    LiftMap(const LiftMap &rhs)
        : lmap(rhs.lmap), s1_map(rhs.s1_map), s2_map(rhs.s2_map) {}

    // move constructor
    LiftMap(LiftMap &&rhs)
        : lmap(std::move(rhs.lmap)),
          s1_map(std::move(rhs.s1_map)),
          s2_map(std::move(s2_map)) {}

    // copy assignment operator
    LiftMap &operator=(const LiftMap &rhs) {
        lmap = rhs.lmap;
        s1_map.clear();
        s2_map.clear();
        return *this;
    }

    // move assignment operator
    LiftMap &operator=(LiftMap &&rhs) {
        lmap = std::move(rhs.lmap);
        s1_map = std::move(s1_map);
        s2_map = std::move(s2_map);
        return *this;
    }

    /* creates a liftover from specified sample in VCF file.
     * input: vcfFile, bcf_hdr_t* of input VCF (via htslib)
     *        sample name of desired individual within the VCF
     * assumes that the vcfFile is sorted!
     */
    LiftMap(vcfFile *fp, bcf_hdr_t *hdr, std::string sample_name,
            std::string haplotype,
            std::vector<std::pair<std::string, std::string>> nm = {},
            LengthMap ls = {}) {
        // std::unordered_map<std::string,size_t> ls = {}) {
        bool get_names = 1;
        if (nm.size()) {
            for (int i = 0; i < nm.size(); ++i)
                add_names(nm[i].first, nm[i].second, i);
            get_names = 0;
        }
        int lmap_idx = 0;  // current index of lmap

        sdsl::bit_vector ibv, dbv, sbv;  // ins, del, snp
        bool no_sample = (sample_name == "");
        if (no_sample)
            fprintf(stderr,
                    "no sample given, assuming GT=1 for all variants\n");
        if (!no_sample && bcf_hdr_set_samples(hdr, sample_name.c_str(), 0))
            die("error: sample does not exist!\n");
        bcf1_t *rec = bcf_init();
        size_t x = 0;
        int rid = -1;
        int ppos = 0;
        int tppos = 0;
        size_t l = 0;
        int hap = stoi(haplotype);
        while (!bcf_read(fp, hdr, rec)) {
            bcf_unpack(rec, BCF_UN_FMT);
            if (rec->rid != rid) {
                if (rid != -1) {  // save the bitvector
                    // condense the bit vectors and add to map
                    if (ppos > l) {
                        die("something went wrong, we went past bounds of s1 "
                            "sequence. "
                            "exiting...\n");
                    }
                    ibv.resize(x + (l - ppos));
                    dbv.resize(x + (l - ppos));
                    sbv.resize(x + (l - ppos));
                    lmap.push_back(Lift(ibv, dbv, sbv));
                }
                rid = rec->rid;
                if (get_names) {
                    const char *name = bcf_hdr_id2name(hdr, rid);
                    add_names(name, name, lmap_idx++);
                }
                // get size of contig
                l = hdr->id[BCF_DT_CTG][rid].val->info[0];
                if (l == 0) {
                    if (!ls.size()) {
                        die("contig length not specified in vcf, nor was it "
                            "provided by "
                            "user!\n");
                    } else {
                        /* get length from input lengths */
                        // auto lpair = ls.find(name_map[bcf_hdr_id2name(hdr,
                        // rid)]); if (lpair != ls.end()) l = lpair->second;
                        // else die("contig length not given by user!!\n");
                        l = ls[rid].second;
                    }
                }
                // initialize bit vectors for this contig
                // TODO: maybe set vec size to something less than 2l
                ibv = sdsl::bit_vector(l * 2);
                dbv = sdsl::bit_vector(l * 2);
                sbv = sdsl::bit_vector(l * 2);
                x = 0;
                ppos = 0;
                tppos = 0;  // end of last variant processed
            }
            int32_t *gt_arr = NULL;
            int32_t ngt_arr = 0;
            int ngt;
            ngt =
                no_sample ? 0 : bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
            int prev_is_ins = 0;
            // only process the variant if it's in the sample's genotype
            if (no_sample || (ngt > 0 && !bcf_gt_is_missing(gt_arr[hap]) &&
                              bcf_gt_allele(gt_arr[hap]))) {
                int var = no_sample ? 1 : bcf_gt_allele(gt_arr[hap]);
                int var_type = bcf_get_variant_type(rec, var);
                int rlen = strlen(rec->d.allele[0]);
                int alen = strlen(rec->d.allele[var]);

                // Determine if overlap. Logic copied from bcftools consensus`:
                // https://github.com/samtools/bcftools/blob/df43fd4781298e961efc951ba33fc4cdcc165a19/consensus.c#L579

                // For some variant types POS+REF refer to the base *before* the
                // event; in such case set trim_beg
                int trim_beg = 0;
                int var_len = rec->d.var[var].n;
                if (var_type & VCF_INDEL)
                    trim_beg = 1;
                else if ((var_type & VCF_OTHER) &&
                         !strcasecmp(rec->d.allele[var], "<DEL>")) {
                    trim_beg = 1;
                    var_len = 1 - rec->rlen;
                } else if ((var_type & VCF_OTHER) &&
                           !strncasecmp(rec->d.allele[var], "<INS", 4))
                    trim_beg = 1;

                if (rec->pos <= tppos) {
                    int overlap = 0;
                    if (rec->pos < tppos || !trim_beg || var_len == 0 ||
                        prev_is_ins)
                        overlap = 1;
                    if (overlap) {
                        if (verbose)
                            fprintf(stderr, "Skipping variant %s:%lld\n",
                                    bcf_seqname(hdr, rec), rec->pos + 1);
                        continue;
                    }
                }

                // at this point, ppos and x should point to *end* of last
                // variant in s1 & s2 alignment, resp.
                x += (rec->pos - ppos);  // x should now be pointing to *start*
                                         // of current variant wrt alignment
                ppos = rec->pos;  // ppos now pointing to *start* of current
                                  // variant
                if (var_type == VCF_INDEL) {
                    ++ppos;
                    ++x;
                    prev_is_ins = 0;
                    if (rlen < alen) {  // ins
                        for (size_t i = 0; i < alen - rlen; ++i) {
                            ibv[x++] = 1;
                        }
                        prev_is_ins = 1;
                    } else if (rlen > alen) {  // del
                        for (size_t i = 0; i < rlen - alen; ++i) {
                            dbv[x++] = 1;
                            ++ppos;  // s1 advances, so ppos advances
                        }
                    }
                    --x;
                    --ppos;
                } else if (var_type == VCF_SNP) {
                    sbv[x] = 1;  // no need to increment x here?
                }
                tppos = rec->pos + rec->rlen -
                        1;  // tppos points to end of this variant
            } else {
                continue;
            }
        }
        if (ppos > l) {
            die("something went wrong, we went past bounds of s1 sequence. "
                "exiting...\n");
        }
        ibv.resize(x + (l - ppos));
        dbv.resize(x + (l - ppos));
        sbv.resize(x + (l - ppos));
        lmap.push_back(Lift(ibv, dbv, sbv));
        length_map = ls;
    }

    // converts a position in the s2 sequence to corresponding position in s1
    // sequence input: sequence name, position within s2 sequence output:
    // position within s1 sequence if liftover is not defined, then returns the
    // original position
    size_t lift_pos(  // s2_to_s1
        std::string n, size_t i) {
        auto it = s2_map.find(n);
        if (it != s2_map.end()) {
            return lmap[it->second].lift_pos(i);
        } else
            return i;
    }

    std::string lift_contig(std::string n, size_t pos) {  // get_other_name()
        if (name_map.find(n) != name_map.end()) return name_map[n];
        // return original name if mapping is not available
        else
            return n;
    }

    /* Convert a CIGAR string in an s2 alignment to the corresponding CIGAR
     * string in the s1 alignment
     *
     * Input: bam1_t alignment record against s2 sequence
     *
     * The bam1_t struct is updated if there's a change in its CIGAR.
     * bam1_t->core.n_cigar specifies numbers of operators in a CIGAR;
     * bam_get_cigar(bam1_t) points to the occurrences.
     *
     * If liftover is not defined, do not update the bam1_t struct.
     */
    void lift_cigar(const std::string &n, bam1_t *b) {
        auto it = s2_map.find(n);
        if (it != s2_map.end()) {
            lmap[it->second].lift_cigar(b);
        }
        // If not found, no update.
    }

    // Note that a mapped read is always liftable
    // under the VcfMap framework.
    void lift_aln(bam1_t *aln, bam_hdr_t *hdr,
                  bam_hdr_t *hdr_dest,  // unused
                  std::string &dest_contig) {
        bam1_core_t c = aln->core;
        size_t pos;
        std::string source_contig;
        size_t lift_status;

        // Paired
        if (c.flag & BAM_FPAIRED) {
            // R1 unmapped
            if (c.flag & BAM_FUNMAP) {
                if (c.flag & BAM_FMUNMAP) {
                    lift_status = LIFT_R1_UM_R2_UM;
                } else {
                    source_contig = hdr->target_name[c.mtid];
                    dest_contig = this->lift_contig(source_contig, c.mpos);
                    pos = this->lift_pos(source_contig, c.mpos);
                    aln->core.pos = pos;
                    lift_status = LIFT_R1_UM_R2_L;
                }
                // R1 mapped
            } else {
                source_contig = hdr->target_name[c.tid];
                dest_contig = this->lift_contig(source_contig, c.pos);
                pos = this->lift_pos(source_contig, c.pos);
                this->lift_cigar(source_contig, aln);
                aln->core.pos = pos;
                if (c.flag & BAM_FMUNMAP) {
                    aln->core.mpos = aln->core.pos;
                    lift_status = LIFT_R1_L_R2_UM;
                } else {
                    std::string msource_contig(hdr->target_name[c.mtid]);
                    std::string mdest_contig(
                        this->lift_contig(msource_contig, c.mpos));
                    c.mtid = sam_hdr_name2tid(hdr, mdest_contig.c_str());
                    size_t mpos = this->lift_pos(msource_contig, c.mpos);
                    aln->core.mpos = mpos;
                    int isize = (c.isize == 0)
                                    ? 0
                                    : c.isize + (mpos - c.mpos) - (pos - c.pos);
                    aln->core.isize = isize;
                    lift_status = LIFT_R1_L_R2_L;
                }
            }
            // Unpaired
        } else {
            if (c.flag & BAM_FUNMAP) {
                lift_status = LIFT_R_UM;
            } else {
                source_contig = hdr->target_name[c.tid];
                dest_contig = this->lift_contig(source_contig, c.pos);
                pos = this->lift_pos(source_contig, c.pos);
                this->lift_cigar(source_contig, aln);
                aln->core.pos = pos;
                lift_status = LIFT_R_L;
            }
        }
    }

    // saves to stream
    size_t serialize(std::ofstream &out) {
        size_t size = 0;
        size_t nelems = lmap.size();
        out.write(reinterpret_cast<char *>(&nelems), sizeof(nelems));
        size += sizeof(nelems);
        for (auto i = 0; i < nelems; ++i) {
            size += lmap[i].serialize(out);
        }
        size += s1_map.serialize(out);
        size += s2_map.serialize(out);
        size += name_map.serialize(out);
        size += LevioSamUtils::serialize_lengthmap(out, length_map);
        return size;
    }

    // loads from stream
    void load(std::ifstream &in) {
        size_t nelems;
        in.read(reinterpret_cast<char *>(&nelems), sizeof(nelems));
        for (auto i = 0; i < nelems; ++i) {
            lmap.push_back(Lift(in));
        }
        s1_map.load(in);
        s2_map.load(in);
        name_map.load(in);
        length_map = LevioSamUtils::load_lengthmap(in);
    }

    // get names and lengths of s1 sequences
    std::pair<std::vector<std::string>, std::vector<size_t>> get_s1_lens() {
        std::vector<size_t> lengths;
        std::vector<std::string> names;
        for (const auto &x : s1_map) {
            lengths.push_back(lmap[x.second].s1_len());
            names.push_back(x.first);
        }
        return std::make_pair(names, lengths);
    }

    LengthMap &get_lengthmap() { return length_map; }

   private:
    std::vector<Lift> lmap;
    LengthMap length_map;

    Name2NameMap name_map;
    Name2IdxMap s1_map;
    Name2IdxMap s2_map;
    int verbose = 0;
    void add_names(std::string n1, std::string n2, int i) {
        if (s1_map.find(n1) == s1_map.end())
            s1_map[n1] = i;
        else
            die("non one-to-one mappings not supported yet!\n");
        if (s2_map.find(n2) == s2_map.end())
            s2_map[n2] = i;
        else
            die("non one-to-one mappings not supported yet!\n");
        if (name_map.find(n1) == name_map.end())
            name_map[n1] = n2;
        else
            die("non one-to-one mappings not supported yet!\n");
        if (n1 != n2) {
            if (name_map.find(n2) == name_map.end())
                name_map[n2] = n1;
            else
                die("non one-to-one mappings not supported yet! Or there are "
                    "identical "
                    "but non-corresponding names between s1 and s2\n");
        }
    }
};

LiftMap lift_from_vcf(std::string fname, std::string sample,
                      std::string haplotype, NameMap names, LengthMap lengths);

};  // namespace lift

#endif
