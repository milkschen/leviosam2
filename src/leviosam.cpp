/*
 * liftover.cpp
 *
 * Classes and routines for translating (lifting over) coordinates
 * between two aligned sequences
 *
 * Authors: Taher Mun, Nae-Chyun Chen, Ben Langmead
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "leviosam.hpp"

#include <getopt.h>
#include <htslib/kseq.h>
#include <htslib/kstring.h>
#include <stdio.h>
#include <zlib.h>

#include <algorithm>
#include <ctime>
#include <vector>

#include "aln.hpp"
#include "collate.hpp"
#include "ksw2.h"
#include "lift_bed.hpp"
#include "reconcile.hpp"
#include "yaml.hpp"

KSEQ_INIT(gzFile, gzread);

/** Updates allocated threads given input args.
 *
 * LevioSAM2 expects either (1) --threads (-t) or (2) --hts_threads and/or
 * --lift_threads is set to a non-default value. If --threads and any of
 * --hts_threads and --lift_threads is set, exits with an error message.
 *
 * If only --threads is set, dynamically assigns --hts_threads to
 * max(1, --threads/4) and --lift_threads to be --threads - --hts_threads.
 * Or updates --threads to be the sum of --hts_threads and --lift_threads.
 *
 * @param args leviosam2 arguments
 * @return 0 if successful, 1 if not.
 */
uint8_t update_thread_allocation(lift_opts &args) {
    if (args.lift_threads < 1) {
        std::cerr << "[E::update_thread_allocation] --lift_threads must "
                     ">= 1.\n";
        return 1;
    }
    if (args.threads > DEFAULT_NUM_THREADS) {
        if (args.lift_threads != DEFAULT_NUM_LIFT_THREADS ||
            args.hts_threads != DEFAULT_NUM_HTS_THREADS) {
            std::cerr << "[E::update_thread_allocation] There are two ways to "
                         "set the number of threads: (1) set --threads (-t). "
                         "With this approach, leviosam2 assigns `max(1, t/4)` "
                         "threads to (de)compress I/O files and the rest to "
                         "the lift core. (2) set --hts_threads and/or "
                         "--lift_threads separately.\n";
            return 1;
        }
        args.hts_threads = std::max(1, args.threads / 4);
        args.lift_threads = args.threads - args.hts_threads;
    } else if (args.lift_threads != DEFAULT_NUM_LIFT_THREADS ||
               args.hts_threads != DEFAULT_NUM_HTS_THREADS) {
        args.threads = args.hts_threads + args.lift_threads;
    }
    std::cerr << "[I::update_thread_allocation] --theads=" << args.threads
              << ", --lift_threads=" << args.lift_threads
              << ", --hts_threads=" << args.hts_threads << "\n";
    return 0;
}

NameMap parse_name_map(const char *fname) {
    NameMap names;
    FILE *fp = fopen(fname, "r");
    char n1[255];
    char n2[255];
    int x;
    while ((x = fscanf(fp, "%[^\t]\t%[^\n]\n", n1, n2)) != EOF) {
        if (x <= 0) {
            std::cerr << "[E::parse_name_map] Failed to read name map from "
                      << fname << "\n";
            exit(1);
        }
        names.push_back(std::make_pair(std::string(n1), std::string(n2)));
    }
    fclose(fp);
    return names;
}

void serialize_run(lift_opts args) {
    if (args.outpre == "") {
        std::cerr << "[E::serialize_run] No output prefix specified. Please "
                     "set `-p`.\n";
        print_serialize_help_msg();
        exit(1);
    }
    if (args.length_map.size() == 0) {
        std::cerr << "[E::serialize_run] No length map is found. Please "
                     "set `-F`.\n";
        print_serialize_help_msg();
        exit(1);
    }

    // VcfMap
    if (args.vcf_fname != "") {
        std::string fn_index = args.outpre + ".lft";
        lift::LiftMap l(lift::lift_from_vcf(args.vcf_fname, args.sample,
                                            args.haplotype, args.name_map,
                                            args.length_map));
        std::ofstream o(fn_index, std::ios::binary);
        l.serialize(o);
        std::cerr << "[I::serialize_run] levioSAM VcfMap saved to " << fn_index
                  << "\n";
        // ChainMap
    } else if (args.chain_fname != "") {
        std::string fn_index = args.outpre + ".clft";
        chain::ChainMap cfp(args.chain_fname, args.verbose,
                            args.allowed_cigar_changes, args.length_map);
        std::ofstream o(fn_index, std::ios::binary);
        cfp.serialize(o);
        std::cerr << "[I::serialize_run] levioSAM ChainMap saved to "
                  << fn_index << "\n";
    } else {
        std::cerr << "[E::serialize_run] Cannot build a levioSAM index. "
                     "Please set -v or -c properly\n";
        print_serialize_help_msg();
        exit(1);
    }
}

template <class T>
void read_and_lift(T *lift_map, std::mutex *mutex_fread,
                   std::mutex *mutex_fwrite, samFile *sam_fp,
                   samFile *out_sam_fp, sam_hdr_t *hdr_source,
                   sam_hdr_t *hdr_dest,
                   const std::map<std::string, std::string> *ref_dict,
                   LevioSamUtils::WriteDeferred *wd, const lift_opts &args) {
    std::string current_contig;
    std::vector<bam1_t *> aln_vec, aln_vec_clone;
    for (int i = 0; i < args.chunk_size; i++) {
        bam1_t *aln = bam_init1();
        aln_vec.push_back(aln);
        bam1_t *aln_clone = bam_init1();
        aln_vec_clone.push_back(aln_clone);
    }
    int read = 1;
    std::string ref;
    while (read >= 0) {
        int num_actual_reads = args.chunk_size;
        {
            // read from SAM, protected by mutex
            std::lock_guard<std::mutex> g(*mutex_fread);
            for (int i = 0; i < args.chunk_size; i++) {
                read = sam_read1(sam_fp, hdr_source, aln_vec[i]);
                if (read < 0) {
                    num_actual_reads = i;
                    break;
                }
            }
        }
        for (int i = 0; i < num_actual_reads; i++) {
            auto clone = bam_copy1(aln_vec_clone[i], aln_vec[i]);
            // TODO: Plan to remove `null_dest_contig`.
            // Need to test scenarios with `NameMap`
            std::string null_dest_contig;
            lift_map->lift_aln(aln_vec[i], hdr_source, hdr_dest,
                               null_dest_contig, args.keep_mapq);
            // Skip update MD, NM tags and realignment if `md_flag` is not
            // set
            if (!args.md_flag) {
                // strip MD and NM tags if `md_flag` not set bc liftover
                // invalidates them
                LevioSamUtils::remove_nm_md_tag(aln_vec[i]);
                continue;
            }
            // if cannot find a chromosome id, remove tags and skip
            int32_t tid = aln_vec[i]->core.tid;
            if (tid == -1) {
                LevioSamUtils::remove_nm_md_tag(aln_vec[i]);
                continue;
            }
            std::string dest_contig(
                hdr_dest->target_name[aln_vec[i]->core.tid]);
            // change ref if needed
            if (dest_contig != current_contig) {
                ref = ref_dict->at(dest_contig);
                current_contig = dest_contig;
            }
            bam_fillmd1(aln_vec[i], ref.data(), args.md_flag, 1);

            // Tags specific to Ultima Genomics reads
            if (args.ultima_genomics == true) {
                bool reversely_lifted = LevioSamUtils::is_reversedly_lifted(
                    aln_vec_clone[i], aln_vec[i]);
                LevioSamUtils::update_ultima_genomics_tags(aln_vec[i],
                                                           reversely_lifted);
            }

            // Skip re-alignment if `-x` is not set
            if (args.realign_yaml == "") continue;
            uint8_t *nm_ptr = bam_aux_get(aln_vec[i], "NM");
            // Skip is the NM tag is not found
            if (nm_ptr == NULL) continue;
            // Skip if unmapped
            if (aln_vec[i]->core.flag & BAM_FUNMAP) continue;
            // Skip if edit distance not greater than the threshold
            if (bam_aux2i(nm_ptr) <= args.aln_opts.nm_threshold) continue;
            std::string ref_seq =
                ref.substr(aln_vec[i]->core.pos,
                           bam_cigar2rlen(aln_vec[i]->core.n_cigar,
                                          bam_get_cigar(aln_vec[i])));
            std::string q_seq = LevioSamUtils::get_read(aln_vec[i]);
            std::vector<uint32_t> new_cigar;
            // We perform global alignment for local sequences and
            // thus an alignment with a high number of clipped bases
            // might not re-align well
            int new_score = 0;
            if (Aln::align_ksw2(ref_seq.data(), q_seq.data(), args.aln_opts,
                                new_cigar, new_score) > 0) {
                Cigar::update_cigar(aln_vec[i], new_cigar);
                // Redo fillmd after re-align
                bam_fillmd1(aln_vec[i], ref.data(), args.md_flag, 1);
                bam_aux_append(aln_vec[i], "LR", 'i', 4, (uint8_t *)&new_score);
            } else {
                // If re-alignment fails, set this read to unmapped
                bool first_seg =
                    (aln_vec[i]->core.flag & BAM_FPAIRED)
                        ? (aln_vec[i]->core.flag & BAM_FREAD1) ? true : false
                        : true;
                LevioSamUtils::update_flag_unmap(aln_vec[i], first_seg,
                                                 args.keep_mapq);
                if (args.verbose >= VERBOSE_INFO) {
                    std::cerr << "[I::read_and_lift] Zero-length new cigar for "
                              << bam_get_qname(aln_vec[i])
                              << "; flag = " << aln_vec[i]->core.flag
                              << "; NM:i = " << bam_aux2i(nm_ptr)
                              << "; l_qseq = " << aln_vec[i]->core.l_qseq
                              << ". Set to unmapped.\n";
                }
            }
        }
        for (int i = 0; i < num_actual_reads; i++) {
            // If a read is committed
            if (wd->commit_aln_dest(aln_vec[i]) == true ||
                wd->commit_aln_source(aln_vec_clone[i]) == true) {
                // write to file, thread corruption protected by lock_guard
                std::lock_guard<std::mutex> g_commited(*mutex_fwrite);
                if (sam_write1(out_sam_fp, hdr_dest, aln_vec[i]) < 0) {
                    std::cerr << "[E::read_and_lift] Failed to write record "
                              << bam_get_qname(aln_vec[i]) << "\n";
                    exit(1);
                }
            }

            // If suppress, only write the original alignment to
            // `unliftable` If defer, write to both `unliftable` and `defer`
            if (wd->commit_aln_dest(aln_vec[i]) == false) {
                std::lock_guard<std::mutex> g_deferred(wd->mutex_fwrite);
                // Write deferred reads
                if (wd->commit_aln_source(aln_vec_clone[i]) == false) {
                    // Optionally realign
                    wd->write_deferred_bam(aln_vec[i], hdr_dest);
                }
                // Write suppressed reads
                wd->write_deferred_bam_orig(aln_vec_clone[i]);
            }
        }
    }
    for (int i = 0; i < args.chunk_size; i++) {
        bam_destroy1(aln_vec[i]);
        bam_destroy1(aln_vec_clone[i]);
    }
    aln_vec.clear();
    aln_vec_clone.clear();
}

/* Load a FASTA file and return a map */
std::map<std::string, std::string> load_fasta(std::string ref_name) {
    std::cerr << "[I::load_fasta] Loading FASTA...";
    std::map<std::string, std::string> fmap;
    gzFile fp = gzopen(ref_name.data(), "r");
    if (!fp) {
        std::cerr << "\n[E::load_fasta] Cannot open file " << ref_name << ".\n";
        exit(1);
    }
    kseq_t *seq;
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        fmap[std::string(seq->name.s)] = std::string(seq->seq.s);
    }
    kseq_destroy(seq);
    gzclose(fp);
    std::cerr << "done\n";
    return fmap;
}

void lift_run(lift_opts args) {
    chain::ChainMap chain_map = [&] {
        if (args.chainmap_fname != "") {
            std::cerr << "[I::lift_run] Loading levioSAM index...";
            std::ifstream in(args.chainmap_fname, std::ios::binary);
            return chain::ChainMap(in, args.verbose,
                                   args.allowed_cigar_changes);
        } else if (args.chain_fname != "") {
            std::cerr << "[I::lift_run] Building levioSAM index...";
            if (args.length_map.size() == 0) {
                std::cerr << "[E::lift_run] No length map is found. Please "
                             "set -F properly.\n";
                print_serialize_help_msg();
                exit(1);
            }
            return chain::ChainMap(args.chain_fname, args.verbose,
                                   args.allowed_cigar_changes, args.length_map);
        } else {
            return chain::ChainMap();
        }
    }();
    lift::LiftMap lift_map = [&] {
        if (args.lift_fname != "") {
            std::cerr << "[I::lift_run] Loading levioSAM index...";
            std::ifstream in(args.lift_fname, std::ios::binary);
            return lift::LiftMap(in);
            // if "-l" not specified, then create a levioSAM
        } else if (args.vcf_fname != "") {
            std::cerr << "[I::lift_run] Building levioSAM index...";
            return lift::LiftMap(
                lift::lift_from_vcf(args.vcf_fname, args.sample, args.haplotype,
                                    args.name_map, args.length_map));
        } else if ((args.chain_fname != "") || (args.chainmap_fname != "")) {
            return lift::LiftMap();
        } else {
            std::cerr << "[E::lift_run] Not enough parameters specified to "
                         "build/load lift-over\n";
            print_lift_help_msg();
            exit(1);
        }
    }();
    if (args.chainmap_fname != "") {
        args.length_map = chain_map.length_map;
    } else if (args.lift_fname != "") {
        args.length_map = lift_map.get_lengthmap();
    }
    if (args.length_map.size() == 0) {
        std::cerr << "[E::lift_run] Length map is empty\n";
        exit(1);
    }
    std::cerr << "done\n";

    samFile *sam_fp = (args.sam_fname == "")
                          ? sam_open("-", "r")
                          : sam_open(args.sam_fname.data(), "r");

    // Gets extra threads to decompress HTS files when requested
    if (args.hts_threads > 0) hts_set_threads(sam_fp, args.hts_threads);

    if (!sam_fp) {
        std::cerr << "[E::lift_run] Invalid alignment input\n";
        exit(1);
    }
    std::string out_mode = (args.out_format == "sam") ? "w" : "wb";
    // if split, append "-committed" after `outpre`
    std::string out_sam_fname =
        (args.split_rules.size() == 0)
            ? args.outpre + "." + args.out_format
            : args.outpre + "-committed." + args.out_format;
    samFile *out_sam_fp = (args.outpre == "-" || args.outpre == "")
                              ? sam_open("-", out_mode.data())
                              : sam_open(out_sam_fname.data(), out_mode.data());
    if (!out_sam_fp) {
        std::cerr << "[E::lift_run] Invalid alignment output\n";
        exit(1);
    }
    // Gets extra threads to compress HTS files when requested
    if (args.hts_threads > 0) hts_set_threads(out_sam_fp, args.hts_threads);

    sam_hdr_t *hdr_orig = sam_hdr_read(sam_fp);
    if (!hdr_orig) {
        std::cerr << "[E::lift_run] Invalid input header\n";
        exit(1);
    }
    sam_hdr_t *hdr = LevioSamUtils::lengthmap_to_hdr(args.length_map, hdr_orig);
    sam_hdr_add_pg(hdr, "leviosam2", "VN", VERSION, "CL", args.cmd.data(),
                   NULL);
    sam_hdr_add_pg(hdr_orig, "leviosam2", "VN", VERSION, "CL", args.cmd.data(),
                   NULL);
    auto write_hdr = sam_hdr_write(out_sam_fp, hdr);

    if (args.md_flag) {
        if (args.ref_name == "") {
            std::cerr << "[E::lift_run] -m/--md -f <fasta> to be "
                         "provided as well\n";
            exit(1);
        }
    }
    std::map<std::string, std::string> ref_dict;
    if (args.ref_name != "") {
        ref_dict = load_fasta(args.ref_name);
        if (ref_dict.size() == 0) {
            std::cerr << "[E::lift_run] Failed to load reference from "
                      << args.ref_name << "\n";
            exit(1);
        }
    }

    LevioSamUtils::WriteDeferred wd;
    if (args.split_rules.size() != 0) {
        wd.init(args.split_rules, hdr_orig, hdr, args.bed_defer_source,
                args.bed_defer_dest, args.bed_commit_source,
                args.bed_commit_dest, args.bed_isec_threshold);
        wd.open_sams(args.outpre, args.out_format);
        wd.print_info();
    }

    std::vector<std::thread> threads;
    std::mutex mutex_fread, mutex_fwrite;
    for (int j = 0; j < args.lift_threads; j++) {
        if (args.chain_fname == "" && args.chainmap_fname == "") {
            threads.push_back(std::thread(read_and_lift<lift::LiftMap>,
                                          &lift_map, &mutex_fread,
                                          &mutex_fwrite, sam_fp, out_sam_fp,
                                          hdr_orig, hdr, &ref_dict, &wd, args));
        } else {
            threads.push_back(std::thread(read_and_lift<chain::ChainMap>,
                                          &chain_map, &mutex_fread,
                                          &mutex_fwrite, sam_fp, out_sam_fp,
                                          hdr_orig, hdr, &ref_dict, &wd, args));
        }
    }
    for (int j = 0; j < args.lift_threads; j++) {
        if (threads[j].joinable()) threads[j].join();
    }
    threads.clear();
    sam_hdr_destroy(hdr);
    sam_hdr_destroy(hdr_orig);
    sam_close(sam_fp);
    sam_close(out_sam_fp);
    wd.close_sams();
}

/* This is a wrapper for the constructor of LiftMap.*/
lift::LiftMap lift::lift_from_vcf(std::string fname, std::string sample,
                                  std::string haplotype, NameMap names,
                                  LengthMap lengths) {
    if (fname == "" && sample == "") {
        std::cerr << "[E::lift_from_vcf] Vcf file name and sample name are "
                     "required.\n";
        print_serialize_help_msg();
        exit(1);
    }
    vcfFile *fp = bcf_open(fname.data(), "r");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    return lift::LiftMap(fp, hdr, sample, haplotype, names, lengths);
}

void print_serialize_help_msg() {
    std::cerr << "\n";
    std::cerr << "Build a levioSAM2 index of a chain file.\n";
    std::cerr << "Usage:   leviosam2 index -c <chain> -p <out_prefix> "
                 "-F <fai>\n";
    std::cerr << "\n";
    std::cerr << "Inputs:  -c path   Path to the chain file to index.\n";
    std::cerr << "         -F path   Path to the FAI (FASTA index) file of the "
                 "target reference.\n";
    std::cerr << "         -p string Prefix of the output file.\n";
    std::cerr << "\n";
}

void print_lift_help_msg() {
    std::cerr << "\n";
    std::cerr << "Lift over using leviosam2.\n";
    std::cerr << "Usage:   leviosam2 lift [options] -C <clft>\n";
    std::cerr << "\n";
    std::cerr << "Inputs:  -C path   Path to an indexed ChainMap.\n";
    std::cerr << "Options: -a path   Path to the SAM/BAM/CRAM file "
                 "to be lifted. \n";
    std::cerr << "                   "
                 "Leave empty or set to \"-\" to read from stdin.\n";
    std::cerr << "         -t INT    Number of threads used.\n"
                 "                   "
                 "If -t is not set, the value would be the sum of\n"
                 "                   "
                 "--hts_threads and --lift_threads. [1] \n";
    std::cerr << "         --lift_threads INT "
                 "Number of threads used for lifting reads. \n"
                 "                            "  // align
                 "If -t is set, the value should be left unset.\n"
                 "                            "  // align
                 "The value would be inferred as `t - max(1, t/4)`. [1]\n";
    std::cerr << "         --hts_threads INT "
                 "Number of threads used to compress/decompress HTS files.\n"
                 "                           "  // align
                 "This can improve thread scaling.\n"
                 "                           "  // align
                 "If -t is set, the value should be left unset.\n"
                 "                           "  // align
                 "The value would be inferred as `max(1, t/4)`. [0]\n";
    std::cerr << "         -m        "
                 "Add MD and NM to output alignment records (requires -f)\n";
    std::cerr << "         -f path   "
                 "Path to the FASTA file of the target reference. \n";
    std::cerr << "         -x path   Re-alignment preset. [] \n";
    std::cerr << "         -G INT    "
                 "Maximum allowed gap size (in base pairs) between chain "
                 "intervals for an alignment to be considered liftable. A "
                 "value of 0 requires perfect chain interval continuity. [0]\n";
    std::cerr << "         -T INT    "
                 "Chunk size for each thread. [256] \n"
                 "                   "  // align
                 "Each thread queries <-T> reads, lifts, and writes.\n"
                 "                   "  // align
                 "Setting a larger -T uses slightly more "
                 "memory but might benefit thread scaling.\n";
    std::cerr << "\n";
    std::cerr << "         Commit/defer rule options:\n";
    std::cerr << "           -S string<:int/float> Key-value pair of "
                 "a split rule. We allow appending multiple `-S` options.\n";
    std::cerr << "                     Options: "
                 "mapq:<int>, aln_score:<int>, isize:<int>, hdist:<int>, "
                 "clipped_frac:<float>. lifted. [none]\n";
    std::cerr << "                       * mapq          INT   "
                 "Min MAPQ (pre-liftover) accepted for a committed read.\n";
    std::cerr << "                       * aln_score     INT   "
                 "Min AS:i (alignment score) (pre-liftover) accepted for a "
                 "committed read.\n";
    std::cerr << "                       * isize         INT   "
                 "Max TLEN/isize (post-liftover) accepted for a committed "
                 "read.\n";
    std::cerr << "                       * hdist         INT   "
                 "Max NM:i (Edit distance) (post-liftover) accepted for a "
                 "committed read.\n"
                 "                                             "
                 "`-m` and `-f` must be set.\n";
    std::cerr << "                       * clipped_frac  FLOAT "
                 "Max fraction of clipped bases (post-liftover) accepted for a "
                 "committed read.\n"
                 "                                             "
                 "A higher value results in fewer committed reads.\n";
    std::cerr << "           Example: `-S mapq:20 -S aln_score:20` commits "
                 "MQ>=20 and AS>=20 alignments.\n";
    std::cerr << "           -r string Path to a BED file (source "
                 "coordinates). Reads overlap with the regions are always "
                 "committed. [none]\n";
    std::cerr << "           -D string Path to a BED file (dest coordinates). "
                 "Reads overlap with the regions are always deferred. [none]\n";
    std::cerr << "           -B float  Threshold for BED record intersection. "
                 "[0]\n"
                 "                     If <= 0: consider any overlap (>0 bp)\n"
                 "                     If > 1: consider >`-B`-bp overlap\n"
                 "                     If 1>=`-B`>0: consider overlap with a "
                 "fraction of >`-B` of the alignment.\n";
    std::cerr << "\n";
}

void print_main_help_msg() {
    std::cerr << "\n";
    std::cerr << "Program: leviosam2 (lift over alignments using a "
                 "chain file)\n";
    std::cerr << "Version: " << VERSION << "\n";
    std::cerr << "Usage:   leviosam2 <command> [options]\n\n";
    std::cerr << "Commands: index       Build a levioSAM2 index of "
                 "a chain file.\n";
    std::cerr << "          lift        Lift alignments in SAM/BAM/CRAM "
                 "formats.\n";
    std::cerr << "          bed         Lift intervals in BED format.\n";
    std::cerr << "          collate     Collate lifted paired-end "
                 "alignments.\n";
    std::cerr << "          reconcile   Reconcile alignments.\n";
    std::cerr << "Options:  -h          Print detailed usage.\n";
    std::cerr << "          -V          Verbose level [0].\n";
    std::cerr << "          --version   Print version.\n";
    std::cerr << "\n";
}

int main(int argc, char **argv) {
    if (argc - optind < 1) {
        std::cerr << "[E::main] No argument provided\n";
        print_main_help_msg();
        exit(1);
    }
    // These functionalities are separated and use different arguments
    if (!strcmp(argv[optind], "collate")) {
        return collate_run(argc, argv);
    } else if (!strcmp(argv[optind], "bed")) {
        return lift_bed_run(argc, argv);
    } else if (!strcmp(argv[optind], "reconcile")) {
        return reconcile_run(argc, argv);
    }

    double start_cputime = std::clock();
    auto start_walltime = std::chrono::system_clock::now();

    int c;
    lift_opts args;
    args.cmd = LevioSamUtils::make_cmd(argc, argv);
    static struct option long_options[]{
        {"sam", required_argument, 0, 'a'},
        {"bed_isec_threshold", required_argument, 0, 'B'},
        {"chain", required_argument, 0, 'c'},
        {"chainmap", required_argument, 0, 'C'},
        {"bed_defer_source", required_argument, 0, 'd'},
        {"bed_defer_dest", required_argument, 0, 'D'},
        {"reference", required_argument, 0, 'f'},
        {"dest_fai", required_argument, 0, 'F'},
        {"haplotype", required_argument, 0, 'g'},
        {"allowed_cigar_changes", required_argument, 0, 'G'},
        {"leviosam", required_argument, 0, 'l'},
        {"md", no_argument, 0, 'm'},
        {"namemap", required_argument, 0, 'n'},
        {"out_format", required_argument, 0, 'O'},
        {"prefix", required_argument, 0, 'p'},
        {"bed_commit_source", required_argument, 0, 'r'},
        {"bed_commit_dest", required_argument, 0, 'R'},
        {"sample", required_argument, 0, 's'},
        {"split_mode", required_argument, 0, 'S'},
        {"threads", required_argument, 0, 't'},
        {"chunk_size", required_argument, 0, 'T'},
        {"vcf", required_argument, 0, 'v'},
        {"verbose", required_argument, 0, 'V'},
        {"realign_yaml", required_argument, 0, 'x'},
        {"keep_mapq", no_argument, 0, OPT_KEEP_MAPQ},
        {"lift_threads", required_argument, 0, OPT_LIFT_THREADS},
        {"hts_threads", required_argument, 0, OPT_HTS_THREADS},
        {"ultima_genomics", no_argument, 0, OPT_ULTIMA_GENOMICS},
        {"help", no_argument, 0, OPT_HELP},
        {"version", no_argument, 0, OPT_VERSION}};
    int long_index = 0;
    while ((c = getopt_long(argc, argv,
                            "hma:B:c:C:d:D:f:F:g:G:l:n:O:p:r:R:s:S:t:T:v:V:x:",
                            long_options, &long_index)) != -1) {
        switch (c) {
            case 'h':
            case OPT_HELP:
                print_main_help_msg();
                print_serialize_help_msg();
                print_lift_help_msg();
                exit(0);
            case 'm':
                args.md_flag |= 8;
                args.md_flag |= 16;
                break;
            case 'a':
                args.sam_fname = optarg;
                break;
            case 'B':
                args.bed_isec_threshold = std::stof(optarg);
            case 'c':
                args.chain_fname = optarg;
                break;
            case 'C':
                args.chainmap_fname = optarg;
                break;
            case 'd':
                std::cerr << "[W::main] -d has not been fully tested\n";
                args.bed_defer_source.init(optarg);
                break;
            case 'D':
                args.bed_defer_dest.init(optarg);
                break;
            case 'f':
                args.ref_name = optarg;
                break;
            case 'F':
                args.length_map = LevioSamUtils::fai_to_map(optarg);
                break;
            case 'g':
                args.haplotype = optarg;
                break;
            case 'G':
                args.allowed_cigar_changes = atoi(optarg);
                break;
            case 'l':
                args.lift_fname = optarg;
                break;
            case 'n':
                args.name_map = parse_name_map(optarg);
                break;
            case 'O':
                args.out_format = optarg;
                break;
            case 'p':
                args.outpre = optarg;
                break;
            case 'r':
                args.bed_commit_source.init(optarg);
                break;
            case 'R':
                std::cerr << "[W::main] -R has not been fully tested\n";
                args.bed_commit_dest.init(optarg);
                break;
            case 's':
                args.sample = optarg;
                break;
            case 'S':
                if (LevioSamUtils::add_split_rule(args.split_rules, optarg) ==
                    false) {
                    exit(1);
                }
                break;
            case 't':
                args.threads = atoi(optarg);
                break;
            case 'T':
                args.chunk_size = atoi(optarg);
                break;
            case 'V':
                args.verbose = atoi(optarg);
                break;
            case 'v':
                args.vcf_fname = optarg;
                break;
            case 'x':
                args.realign_yaml = optarg;
                break;
            case OPT_ULTIMA_GENOMICS:
                args.ultima_genomics = true;
                break;
            case OPT_VERSION:
                std::cout << VERSION << std::endl;
                exit(0);
            case OPT_KEEP_MAPQ:
                args.keep_mapq = true;
                break;
            case OPT_LIFT_THREADS:
                args.lift_threads = atoi(optarg);
                break;
            case OPT_HTS_THREADS:
                args.hts_threads = atoi(optarg);
                break;
            default:
                std::cerr << "[E::main] Invalid flag value " << c << "\n";
                exit(1);
        }
    }

    if (args.haplotype != "0" && args.haplotype != "1") {
        std::cerr << "[E::main] Invalid haplotype " << args.haplotype.c_str()
                  << "\n";
        exit(1);
    }
    if (args.out_format != "sam" && args.out_format != "bam") {
        std::cerr << "[E::main] Not supported extension format "
                  << args.out_format.c_str() << "\n";
        exit(1);
    }

    for (auto &r : args.split_rules) {
        if (r.first == "hdist") {
            if (args.ref_name == "") {
                std::cerr << "[E::main] Option `-f` must be set when `-S "
                             "hdist` is used\n";
                print_lift_help_msg();
                exit(1);
            }
        }
    }
    if (args.split_rules.size() == 0) {
        if (args.bed_defer_source.get_fn() != "" ||
            args.bed_defer_dest.get_fn() != "" ||
            args.bed_commit_source.get_fn() != "" ||
            args.bed_commit_dest.get_fn() != "") {
            std::cerr << "[E::main] `-S` should be set if any among `-d/D/r/R` "
                         "is set.\n";
            exit(1);
        }
    }
    if (args.realign_yaml != "") {
        if (args.ref_name == "") {
            std::cerr << "[E::main] Option `-f` must be non-empty when `-x` "
                         "is not empty\n";
            exit(1);
        } else {
            const char *fc = args.realign_yaml.c_str();
            std::string yaml = Yaml::file_get_contents<std::string>(fc);
            ryml::Tree realign_tree =
                ryml::parse_in_arena(ryml::to_csubstr(yaml));
            args.aln_opts.deserialize_realn(realign_tree);
            args.aln_opts.print_parameters();
        }
    }

    if (update_thread_allocation(args)) {
        exit(ERROR_THREADS_ALLOCATION);
    }

    if (!strcmp(argv[optind], "lift")) {
        lift_run(args);
    } else if (!strcmp(argv[optind], "serialize") ||
               !strcmp(argv[optind], "index")) {
        serialize_run(args);
    }

    double cpu_duration =
        (std::clock() - start_cputime) / (double)CLOCKS_PER_SEC;
    std::chrono::duration<double> wall_duration =
        (std::chrono::system_clock::now() - start_walltime);
    std::cerr << "\n";
    std::cerr << "[I::main] Finished in " << cpu_duration << " CPU seconds, or "
              << wall_duration.count() << " wall clock seconds\n";
    return 0;
}
