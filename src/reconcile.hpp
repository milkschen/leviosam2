#ifndef RECONCILE_H__
#define RECONCILE_H__
#include <tuple>
#include <vector>

#include <htslib/sam.h>


const int MERGE_PE_SUM = 0;
const int MERGE_PE_MAX = 1;
// const int MERGE_PE_PAIR_TLEN_THRESHOLD = 2000;

struct reconcile_opts{
    bool paired_end = false;
    std::string cmd = "";
    int rand_seed = 0;
    std::string output_fn = "";
    std::string score_tag = "NM";
    int merge_pe_mode = MERGE_PE_SUM;
    std::vector<std::string> inputs;
};


template <typename A, typename B, typename C, typename D>
std::vector<std::tuple<A, B, C, D>> zip(const std::vector<A> &a,
                                        const std::vector<B> &b,
                                        const std::vector<C> &c,
                                        const std::vector<D> &d);

template <typename A, typename B, typename C, typename D>
void unzip(const std::vector<std::tuple<A, B, C, D>> &zipped,
           std::vector<A> &a,
           std::vector<B> &b,
           std::vector<C> &c,
           std::vector<D> &d);

/* Random generator function
 *
 * From http://www.cplusplus.com/reference/algorithm/random_shuffle/
 */
static int myrandom (int i) { return std::rand() % i;}

std::vector<std::string> read_file_as_vector(std::string list_fn);

int select_best_aln(const std::vector<bool>& pair_indicators,
                    const std::vector<int>& scores,
                    const std::vector<int>& mapqs,
                    int& num_tied_best);

int select_best_aln_single_end(
    const std::vector<bam1_t*>& aln1s,
    const std::string& score_tag);

int select_best_aln_paired_end(
    const std::vector<bam1_t*>& aln1s,
    const std::vector<bam1_t*>& aln2s,
    const std::string& score_tag,
    const int merge_pe_mode);

void reconcile_core(
    const std::vector<std::string>& sam_fns,
    const std::vector<std::string>& ids,
    const std::vector<samFile*>& sam_fps,
    const std::vector<bam_hdr_t*>& hdrs,
    const bool is_paired_end,
    samFile* out_fp,
    const std::string& score_tag
);

void reconcile(reconcile_opts args);

int reconcile_run(int argc, char** argv);

static void print_reconcile_help();

#endif /* RECONCILE_H__ */
