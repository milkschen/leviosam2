#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <regex>

#include <getopt.h>
#include <limits.h>
#include <string.h>

#include "reconcile.hpp"
#include "leviosam_utils.hpp"

/* Fill the zipped vector with pairs consisting of the corresponding elements of
 * a, b, c and d. (This assumes that the vectors have equal length)
 *
 * Adapted from Marco13:
 * https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one/46370189
 */
template <typename A, typename B, typename C, typename D>
std::vector<std::tuple<A, B, C, D>> zip(const std::vector<A> &a,
                                        const std::vector<B> &b,
                                        const std::vector<C> &c,
                                        const std::vector<D> &d) {
    std::vector<std::tuple<A, B, C, D>> zipped;
    for (int i = 0; i < a.size(); ++i) {
        zipped.push_back(std::make_tuple(a[i], b[i], c[i], d[i]));
    }
    return zipped;
}


/* Write the first and second element of the pairs in the given zipped vector into
 * a, b, c and d. (This assumes that the vectors have equal length)
 *
 * Adapted from Marco13:
 * https://stackoverflow.com/questions/37368787/c-sort-one-vector-based-on-another-one/46370189
 */
template <typename A, typename B, typename C, typename D>
void unzip(const std::vector<std::tuple<A, B, C, D>> &zipped,
           std::vector<A> &a,
           std::vector<B> &b,
           std::vector<C> &c,
           std::vector<D> &d) {
    for (int i = 0; i < a.size(); i++) {
        a[i] = std::get<0>(zipped[i]);
        b[i] = std::get<1>(zipped[i]);
        c[i] = std::get<2>(zipped[i]);
        d[i] = std::get<3>(zipped[i]);
    }
}


/* Rank a set of alignments.
 * Order by: is_proper_pair > score > MAPQ
 */
int select_best_aln(const std::vector<bool>& pair_indicators,
                    const std::vector<int>& scores,
                    const std::vector<int>& mapqs,
                    int& num_tied_best) {
    int vec_size = pair_indicators.size();
    std::vector<int> ranks(vec_size);
    std::iota(ranks.begin(), ranks.end(), 0);

    std::vector<std::tuple<int, int, bool, int>> zipped = zip(mapqs, scores, pair_indicators, ranks);
    std::sort(zipped.begin(), zipped.end(),
        [&](const auto& a, const auto& b)
        {
            return (std::get<2>(a) != std::get<2>(b))? std::get<2>(a) :
                   (std::get<1>(a) != std::get<1>(b))? (std::get<1>(a) > std::get<1>(b)) :
                   (std::get<0>(a) > std::get<0>(b));
        });
    num_tied_best = vec_size;
    for (int i = 0; i < vec_size - 1; i++) {
        // Compare elements in tuples (ranks should be excluded so cannot simply compare tuples).
        if (std::get<2>(zipped[i]) != std::get<2>(zipped[i+1]) ||
            std::get<1>(zipped[i]) != std::get<1>(zipped[i+1]) ||
            std::get<0>(zipped[i]) != std::get<0>(zipped[i+1])) {
            num_tied_best = i + 1;
            break;
        }
    }
    ranks.clear();

    for (int i = 0; i < vec_size; i++) {
        ranks.push_back(std::get<3>(zipped[i]));
    }
    std::random_shuffle(ranks.begin(), ranks.begin() + num_tied_best, myrandom);

    return ranks[0];
}


int select_best_aln_single_end(
    const std::vector<bam1_t*>& aln1s,
    const std::string& score_tag
) {
    std::vector<int> mapqs, scores;
    // We don't actually need `pair_indicators` in single-end mode.
    // Create this to make it easier to re-use the core comparison function.
    std::vector<bool> pair_indicators;
    for (int i = 0; i < aln1s.size(); i++) {
        bam1_core_t c_aln1 = aln1s[i]->core;
        pair_indicators.push_back(true);
        mapqs.push_back(c_aln1.qual);
        uint8_t* tag1 = bam_aux_get(aln1s[i], score_tag.data());
        int score = ((c_aln1.flag & BAM_FUNMAP) || (tag1 == NULL))?
            INT_MIN : -bam_aux2i(tag1);
        if (!(c_aln1.flag & BAM_FUNMAP) && (tag1 == NULL)) {
            std::cerr << "[W::select_best_aln_paired_end] "
                      << score_tag << " tag is missing for read "
                      << bam_get_qname(aln1s[i]) << "\n";
        }
        // int score = (c_aln1.flag & BAM_FUNMAP)? INT_MIN :
        //                                         -bam_aux2i(bam_aux_get(aln1s[i], "NM"));
        scores.push_back(score);
    }
    int num_tied_best;
    int best_idx = select_best_aln(
        pair_indicators=pair_indicators, scores=scores,
        mapqs=mapqs, num_tied_best=num_tied_best);

    return best_idx;
}


int select_best_aln_paired_end(
    const std::vector<bam1_t*>& aln1s,
    const std::vector<bam1_t*>& aln2s,
    const std::string& score_tag,
    const int merge_pe_mode = MERGE_PE_SUM
) {
    std::vector<int> mapqs, scores;
    // Indicator 
    //      1 if a pair is properly paired (TLEN != 0). We apply a loose threshold here.
    //      0 otherwise
    std::vector<bool> pair_indicators;
    for (int i = 0; i < aln1s.size(); i++) {
        bam1_core_t c_aln1 = aln1s[i]->core, c_aln2 = aln2s[i]->core;
        pair_indicators.push_back(c_aln1.isize != 0);
        if (merge_pe_mode == MERGE_PE_SUM) {
            // MERGE_PE_SUM mode sums MAPQ and AS. AS is set to 0 for an unaligned read.
            mapqs.push_back(c_aln1.qual + c_aln2.qual);
            int score = 0;
            uint8_t* tag1 = bam_aux_get(aln1s[i], score_tag.data());
            int score1 = INT_MIN / 2;
            score += ((c_aln1.flag & BAM_FUNMAP) || (tag1 == NULL))?
                INT_MIN / 2 : -bam_aux2i(tag1);
            if (!(c_aln1.flag & BAM_FUNMAP) && (tag1 == NULL)) {
                std::cerr << "[W::select_best_aln_paired_end] NM tag is missing for read " << bam_get_qname(aln1s[i]) << "\n";
            }

            uint8_t* tag2 = bam_aux_get(aln2s[i], score_tag.data());
            score += ((c_aln2.flag & BAM_FUNMAP) || (tag2 == NULL))?
                INT_MIN / 2 : -bam_aux2i(tag2);
            if (!(c_aln2.flag & BAM_FUNMAP) && (tag2 == NULL)) {
                std::cerr << "[W::select_best_aln_paired_end] NM tag is missing for read " << bam_get_qname(aln2s[i]) << "\n";
            }
            scores.push_back(score);
        } else if (merge_pe_mode == MERGE_PE_MAX) {
            // MERGE_PE_MAX mode takes max MAPQ and AS.
            if (c_aln1.qual > c_aln2.qual)
                mapqs.push_back(c_aln1.qual);
            else
                mapqs.push_back(c_aln2.qual);
            int score = INT_MIN;
            if (!(c_aln1.flag & BAM_FUNMAP))
                score = -bam_aux2i(bam_aux_get(aln1s[i], score_tag.data()));
            if (!(c_aln2.flag & BAM_FUNMAP))
                score = (score > -bam_aux2i(bam_aux_get(aln2s[i], score_tag.data())))?
                    score : -bam_aux2i(bam_aux_get(aln2s[i], score_tag.data()));
            scores.push_back(score);
        } else{
            std::cerr << "[E::select_best_aln_paired_end] Invalid merging mode for paired-end alignments "
                      << merge_pe_mode << "\n";
        }
    }
    int num_tied_best;
    int best_idx = select_best_aln(
        pair_indicators=pair_indicators, scores=scores, mapqs=mapqs, num_tied_best=num_tied_best);

    return best_idx;
}


void reconcile_core(
    const std::vector<std::string>& sam_fns,
    const std::vector<std::string>& ids,
    const std::vector<samFile*>& sam_fps,
    const std::vector<bam_hdr_t*>& hdrs,
    const bool is_paired_end,
    samFile* out_fp,
    const std::string& score_tag
) {
    std::vector<bam1_t*> aln1s, aln2s;
    for (int i = 0; i < sam_fns.size(); i++) {
        aln1s.push_back(bam_init1());
        if (is_paired_end)
            aln2s.push_back(bam_init1());
    }
    bool end = false;
    int num_records = 0;
    while (!end) {
        // If in paired-end mode: read two reads from each of the SAM files in each iteration.
        for (int i = 0; i < sam_fns.size(); i++) {
            if (is_paired_end) {
                while (1) {
                    int read1 = sam_read1(sam_fps[i], hdrs[i], aln1s[i]);
                    if (read1 < 0) {
                        end = true;
                        break;
                    }
                    if (!(aln1s[i]->core.flag & BAM_FSECONDARY) && 
                        !(aln1s[i]->core.flag & BAM_FSUPPLEMENTARY)) {
                        break;
                    }
                }
                while (1) {
                    int read2 = sam_read1(sam_fps[i], hdrs[i], aln2s[i]);
                    if (read2 < 0) {
                        end = true;
                        break;
                    }
                    if (!(aln2s[i]->core.flag & BAM_FSECONDARY) &&
                        !(aln2s[i]->core.flag & BAM_FSUPPLEMENTARY)) {
                        break;
                    }
                }
                // Check read names: they should be identical.
                if (strcmp(bam_get_qname(aln1s[i]), bam_get_qname(aln2s[i])) != 0) {
                    std::cerr << "[E::reconcile_core] SAM file should be sorted by read name.\n";
                    std::cerr << "This can be done using `samtools sort -n`\n";
                    std::cerr << "Mismatched reads: " << bam_get_qname(aln1s[i]) <<
                                 " and " << bam_get_qname(aln2s[i]) << "\n";
                    exit(1);
                }
                if (end) {
                    break;
                }
            } else {
                // Single-end mode
                if (sam_read1(sam_fps[i], hdrs[i], aln1s[i]) < 0) {
                    end = true;
                    break;
                }
            }
            num_records ++;

            // Check if read names from all files are identical.
            if (i > 0)
                if (strcmp(bam_get_qname(aln1s[0]), bam_get_qname(aln1s[i])) != 0) {
                    std::cerr << "[E::reconcile_core] Reads mismatch across SAM files.\n";
                    std::cerr << "Mismatched reads: " << bam_get_qname(aln1s[0]) <<
                                 " and " << bam_get_qname(aln1s[i]) << "\n";
                    exit(1);
                }
        }

        if (end)
            break;
        if (is_paired_end) {
            int best_idx = select_best_aln_paired_end(aln1s, aln2s, score_tag, MERGE_PE_SUM);
            bam_aux_append(
                aln1s[best_idx], "RF", 'Z', ids[best_idx].length() + 1,
                reinterpret_cast<uint8_t*>(const_cast<char*>(ids[best_idx].c_str())));
            bam_aux_append(
                aln2s[best_idx], "RF", 'Z', ids[best_idx].length() + 1,
                reinterpret_cast<uint8_t*>(const_cast<char*>(ids[best_idx].c_str())));
            if (sam_write1(out_fp, hdrs[0], aln1s[best_idx]) < 0 ||
                sam_write1(out_fp, hdrs[0], aln2s[best_idx]) < 0) {
                std::cerr << "[E::reconcile_core] Failed to write to file.\n";
                std::cerr << bam_get_qname(aln1s[best_idx]);
                exit(1);
            }
        } else {
            int best_idx = select_best_aln_single_end(aln1s, score_tag);
            bam_aux_append(
                aln1s[best_idx], "RF", 'Z', ids[best_idx].length() + 1,
                reinterpret_cast<uint8_t*>(const_cast<char*>(ids[best_idx].c_str())));
            if (sam_write1(out_fp, hdrs[0], aln1s[best_idx]) < 0) {
                std::cerr << "[E::reconcile_core] Failed to write to file.\n";
                exit(1);
            }
        }
    }
    for (int i = 0; i < sam_fns.size(); i++) {
        bam_destroy1(aln1s[i]);
        if (is_paired_end)
            bam_destroy1(aln2s[i]);
    }

    if (is_paired_end)
        std::cerr << "[I::reconcile_core] Processed " << num_records << " pairs of reads\n";
    else
        std::cerr << "[I::reconcile_core] Processed " << num_records << " reads\n";
}


void reconcile(reconcile_opts args) {
    if (args.paired_end)
        std::cerr << "[I::reconcile] Paired-end mode\n";
    else
        std::cerr << "[I::reconcile] Single-end mode\n";

    std::srand(args.rand_seed);
    std::vector<std::string> sam_fns;
    std::vector<std::string> ids;
    for (auto& s: args.inputs) {
        std::regex regexz(":");
        std::vector<std::string> vec(
            std::sregex_token_iterator(s.begin(), s.end(), regexz, -1),
            std::sregex_token_iterator());
        if (vec.size() != 2) {
            std::cerr << "[E::reconcile] Invalid format: " << s << "\n";
            exit(1);
        }
        ids.push_back(vec[0]);
        sam_fns.push_back(vec[1]);
    }
    for (int i = 0; i < sam_fns.size(); i++) {
        std::cerr << "File " << i << ": ";
        std::cerr << sam_fns[i] << " (" << ids[i] << ")\n";
    }

    std::vector<samFile*> sam_fps;
    std::vector<bam_hdr_t*> hdrs;
    for (int i = 0; i < sam_fns.size(); i++) {
        // Read each SAM file listed in `--sam_list`.
        sam_fps.push_back(sam_open(sam_fns[i].data(), "r"));
        hdrs.push_back(sam_hdr_read(sam_fps[i]));
        if (sam_hdr_nref(hdrs[i]) != sam_hdr_nref(hdrs[0])) {
            std::cerr << "[W::reconcile] Num REF in `" << sam_fns[i] << "` differs with `"
                      << sam_fns[0] << "`. Please check.\n";
        }
    }
    std::string out_fn = args.output_fn;
    std::string out_mode = "wb";
    if (out_fn.substr(out_fn.find_last_of(".") + 1) == "sam") {
        out_mode = "w";
    }
    samFile* out_fp = sam_open(out_fn.data(), out_mode.data());
    if (sam_hdr_write(out_fp, hdrs[0])) {
        std::cerr << "[E::reconcile] Failed to write SAM header to file " << out_fn << "\n";
        exit(1);
    }
    reconcile_core(sam_fns, ids, sam_fps, hdrs,
                     args.paired_end, out_fp, args.score_tag);
    for (auto& s: sam_fps) {
        sam_close(s);
    }
    sam_close(out_fp);

}


static void print_reconcile_help(){
    std::cerr << "\n";
    std::cerr << "Usage: leviosam reconcile [options] -s <label:input> -o <out>\n";
    std::cerr << "\n";
    std::cerr << "  -s <string:string>  Input label and file; separated by a colon, e.g.\n";
    std::cerr << "                      `-s foo:foo.bam -s bar:bar.bam`\n";
    std::cerr << "  -o <string> Path to the output SAM/BAM file\n";
    std::cerr << "\n";
    std::cerr << "Options:\n";
    std::cerr << "  -m          Set to perform merging in pairs [false]\n";
    std::cerr << "  -r <int>    Random seed used by the program [0]\n";
    std::cerr << "\n";
}


int reconcile_run(int argc, char** argv) {
    int c;
    reconcile_opts args;
    args.cmd = LevioSamUtils::make_cmd(argc, argv);
    static struct option long_options[]{
        {"paired_end", no_argument, 0, 'm'},
        {"output_fn", required_argument, 0, 'o'},
        {"rand_seed", required_argument, 0, 'r'},
        {"inputs", required_argument, 0, 's'},
        {"score_tag", required_argument, 0, 'x'}
    };
    int long_idx = 0;
    while ((c = getopt_long(argc, argv, "hmo:r:s:x:", long_options, &long_idx)) != -1) {
        switch (c) {
            case 'm':
                args.paired_end = true;
                break;
            case 'o':
                args.output_fn = optarg;
                break;
            case 'r':
                args.rand_seed = atoi(optarg);
                break;
            case 's':
                args.inputs.push_back(optarg);
                break;
            case 'x':
                args.score_tag = optarg;
                break;
            case 'h':
                print_reconcile_help();
                exit(1);
            default:
                std::cerr << "[E::reconcile_run] Invalid option " << c << " \n";
                print_reconcile_help();
                exit(1);
        }
    }
    if (args.inputs.size() == 0) {
        std::cerr << "[E::reconcile_run] required argument missed.\n";
        print_reconcile_help();
        exit(1);
    }
    reconcile(args);

    return 0;
}

