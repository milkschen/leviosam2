#ifndef LIFTOVER_HPP
#define LIFTOVER_HPP

#include <cstdio>
// #include <memory>
#include <iostream>
#include <htslib/vcf.h>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/util.hpp>
// #include "flat_hash_map.hpp"

namespace lift {

class Lift {

    public:

    Lift() {
        init_rs_sls();
    }

    Lift(std::ifstream& in) {
        this->load(in);
    }

    Lift(sdsl::bit_vector i, sdsl::bit_vector d, sdsl::bit_vector s) : ins(i), del(d), snp(s) {
            // create rank and select data-structures for everything
            init_rs_sls();
    }

    Lift(const Lift& rhs) : ins(rhs.ins), del(rhs.del), snp(rhs.snp) {
        init_rs_sls();
    }

    Lift(Lift&& rhs) : ins(std::move(rhs.ins)), del(std::move(rhs.del)), snp(std::move(rhs.snp)){
        init_rs_sls();
    }

    Lift& operator=(const Lift& rhs) {
        ins = rhs.ins;
        del = rhs.del;
        snp = rhs.snp;
        init_rs_sls();
        return *this;
    }

    Lift& operator=(Lift&& rhs) {
        ins = std::move(rhs.ins);
        del = std::move(rhs.del);
        snp = std::move(rhs.snp);
        init_rs_sls();
        return *this;
    }

    size_t alt_to_ref(size_t p) const {
        auto x = del_sls0(p+1);
        auto y = ins_rs0(x); // counts number of 0s before x, which, coincidentally, gives position in reference
        return y;
        // return ins_rs0(del_sls0(p+1));
    }

    size_t ref_to_alt(size_t p) const {
        return del_rs0(ins_sls0(p+1));
    }

    size_t reflen() {
        return ins[ins.size() - 1] ? ins_rs0(ins.size() - 1) : ins_rs0(ins.size() - 1) + 1;
    }

    size_t altlen() {
        return del[del.size() - 1] ? del_rs0(del.size() - 1) : del_rs0(del.size() - 1) + 1;
    }

    sdsl::sd_vector<> get_ins() const {
        return ins;
    }

    size_t serialize(std::ofstream& out) const {
        size_t size = 0;
        size += ins.serialize(out);
        size += del.serialize(out);
        size += snp.serialize(out);
        return size;
    }

    void load(std::istream& in) {
        ins.load(in);
        del.load(in);
        snp.load(in);
        init_rs_sls();
    }

    private:

    void init_rs_sls() {
        sdsl::util::init_support(ins_rs0, &ins);
        sdsl::util::init_support(del_rs0, &del);
        sdsl::util::init_support(snp_rs0, &snp);
        sdsl::util::init_support(ins_sls0, &ins);
        sdsl::util::init_support(del_sls0, &del);
        sdsl::util::init_support(snp_sls0, &snp);
    }

    sdsl::sd_vector<> ins;
    sdsl::sd_vector<> del;
    sdsl::sd_vector<> snp;
    sdsl::sd_vector<>::rank_0_type ins_rs0;
    sdsl::sd_vector<>::rank_0_type del_rs0;
    sdsl::sd_vector<>::rank_0_type snp_rs0;
    sdsl::sd_vector<>::select_0_type ins_sls0;
    sdsl::sd_vector<>::select_0_type del_sls0;
    sdsl::sd_vector<>::select_0_type snp_sls0;
};

class LiftMap {

    public:

    LiftMap() {}

    /* TODO: save to file */
    LiftMap(std::ifstream& in) {
        this->load(in);
    }

    // copy constructor
    LiftMap(const LiftMap& rhs) : lmap(rhs.lmap), names(rhs.names) {
    }

    // move constructor
    LiftMap(LiftMap&& rhs) : lmap(std::move(rhs.lmap)), names(std::move(rhs.names)) {
    }

    // copy assignment operator
    LiftMap& operator=(const LiftMap& rhs) {
        lmap = rhs.lmap;
        names.clear();
        names = rhs.names;
        return *this;
    }

    // move assignment operator
    LiftMap& operator=(LiftMap&& rhs) {
        lmap = std::move(rhs.lmap);
        names = std::move(rhs.names);
        return *this;
    }

    /* creates a liftover from specified sample in VCF file. 
     * For simplicity, looks at only the first haplotype of the specified
     * sample
     */
    LiftMap(vcfFile* fp, bcf_hdr_t* hdr, std::string sample_name) {
        sdsl::bit_vector ibv, dbv, sbv; // ins, del, snp
        bcf_hdr_set_samples(hdr, sample_name.c_str(), 0);
        bcf1_t* rec = bcf_init();
        size_t x = 0;
        int rid = -1;
        int ppos = 0;
        int tppos = 0;
        size_t l = 0;
        while(!bcf_read(fp, hdr, rec)) {
            bcf_unpack(rec, BCF_UN_FMT);
            if (rec->rid != rid) {
                if (rid != -1) {
                    // condense the bit vectors and add to map
                    if (ppos > l) {
                        fprintf(stderr, "something went wrong, we went past bounds of ref sequence. exiting...\n");
                        exit(1);
                    }
                    ibv.resize(x + (l - ppos));
                    dbv.resize(x + (l - ppos));
                    sbv.resize(x + (l - ppos));
                    // lmap[std::string(bcf_hdr_id2name(hdr, rid))] = Lift(ibv, dbv, sbv);
                    lmap.push_back(Lift(ibv,dbv,sbv));
                    names.push_back(bcf_hdr_id2name(hdr, rid));
                }
                // get size of contig
                rid = rec->rid;
                l = hdr->id[BCF_DT_CTG][rid].val->info[0];
                // initialize bit vectors for this contig
                ibv = sdsl::bit_vector(l*2);
                dbv = sdsl::bit_vector(l*2);
                sbv = sdsl::bit_vector(l*2);
                x = 0;
                ppos = 0;
                tppos = 0; // start of last variant processed
            }
            if (rec->pos < tppos) {
                fprintf(stderr, "VCF not sorted! %s, %d -> %d\n", rec->d.id, tppos, rec->pos);
                exit(1);
            }
            // at this point, ppos and x should point to *end* of last variant in ref & alignment, resp.
            x += (rec->pos - ppos); // x should now be pointing to *start* of current variant wrt alignment
            tppos = rec->pos;
            ppos = rec->pos;  // ppos now pointing to *start* of current variant
            // check if ordered
            int32_t* gt_arr = NULL;
            int32_t ngt_arr = 0;
            int ngt;
            ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &ngt_arr);
            if (ngt > 0 && !bcf_gt_is_missing(gt_arr[0]) && bcf_gt_allele(gt_arr[0])) {
                int var = bcf_gt_allele(gt_arr[0]);
                int var_type = bcf_get_variant_type(rec, var);
                int rlen = strlen(rec->d.allele[0]);
                int alen = strlen(rec->d.allele[var]);
                if (var_type == VCF_INDEL) {
                    ++ppos;
                    ++x; // per VCF format, first base btwn REF and ALT match, rest indicate indel
                    if (rlen < alen) { // ins
                        for (size_t i = 0; i < alen - rlen; ++i) {
                            ibv[x++] = 1;
                        }
                    } else if (rlen > alen) { // del
                        for (size_t i = 0; i < rlen - alen; ++i) {
                            dbv[x++] = 1;
                            ++ppos; // ref advances, so ppos advances
                        }
                    }
                    --x;
                    --ppos;
                } else if (var_type == VCF_SNP) {
                    sbv[x] = 1; // no need to increment x here?
                }
            } else {
                // fprintf(stderr, "skipping non-snp/non-indel variant at %d\n", rec->pos);
            }
        } 
        if (ppos > l) {
            fprintf(stderr, "something went wrong, we went past bounds of ref sequence. exiting...\n");
            exit(1);
        }
        ibv.resize(x + (l - ppos));
        dbv.resize(x + (l - ppos));
        sbv.resize(x + (l - ppos));
        // lmap[bcf_hdr_id2name(hdr, rid)] = Lift(ibv, dbv, sbv); // move assignment operator
        lmap.push_back(Lift(ibv,dbv,sbv));
        names.push_back(bcf_hdr_id2name(hdr, rid));
    }

    size_t alt_to_ref(std::string chr, size_t i) {
        // return lmap[chr].alt_to_ref(i);
        auto it = std::find(names.begin(), names.end(), chr);
        if (it != std::end(names)) {
            size_t index = std::distance(names.begin(), it);
            return lmap[index].alt_to_ref(i);
        } else {
            return -1;
        }
    }

    size_t serialize(std::ofstream& out) {
        size_t size = 0;
        size_t nelems = lmap.size();
        out << nelems;
        for (const auto& it: lmap) {
            // serialize the bit vectors
            size += it.serialize(out);
        }
        std::vector<std::string> strs;
        for (const auto& it: names) {
            // serialize the names
            size += it.size();
            strs.push_back(it);
        }
        for (auto s: strs) {
            out << s << "\t";
        }
        return size;
    }

    void load(std::ifstream& in) {
        int nelems;
        in >> nelems;
        for (auto i = 0; i < nelems; ++i) {
            lmap.push_back(Lift(in));
        }
        for (auto i = 0; i < nelems; ++i) {
            std::string name;
            in >> name;
            names.push_back(name);
        }
    }

    std::pair<std::vector<std::string>, std::vector<size_t>> get_reflens() {
        std::vector<size_t> lengths;
        for (auto i = 0; i < lmap.size(); ++i) {
            lengths.push_back(lmap[i].reflen());
        }
        return std::make_pair(names, lengths);
    }

    std::pair<std::vector<std::string>, std::vector<size_t>> get_altlens() {
        std::vector<size_t> lengths;
        for (auto i = 0; i < lmap.size(); ++i) {
            lengths.push_back(lmap[i].altlen());
        }
        return std::make_pair(names, lengths);
    }

    private:

    // ska::flat_hash_map<std::string, Lift> lmap;
    std::vector<Lift> lmap;
    std::vector<std::string> names;
};

};

#endif
