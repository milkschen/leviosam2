#include "aln.hpp"


namespace Aln {

void AlnOpts::deserialize_realn(ryml::Tree realign_tree) {
    if (realign_tree["alignment"].has_child("engine"))
        realign_tree["alignment"]["engine"] >> engine;
    if (realign_tree["alignment"].has_child("flag"))
        realign_tree["alignment"]["flag"] >> flag;
    if (realign_tree["alignment"].has_child("nm_threshold"))
        realign_tree["alignment"]["nm_threshold"] >> nm_threshold;
    if (realign_tree["alignment"]["a"].has_key())
    if (realign_tree["alignment"].has_child("a"))
        realign_tree["alignment"]["a"] >> a;
    if (realign_tree["alignment"].has_child("b"))
        realign_tree["alignment"]["b"] >> b;
    if (realign_tree["alignment"].has_child("q"))
        realign_tree["alignment"]["q"] >> q;
    if (realign_tree["alignment"].has_child("e"))
        realign_tree["alignment"]["e"] >> e;
    if (realign_tree["alignment"].has_child("q2"))
        realign_tree["alignment"]["q2"] >> q2;
    if (realign_tree["alignment"].has_child("e2"))
        realign_tree["alignment"]["e2"] >> e2;
    if (realign_tree["alignment"].has_child("w"))
        realign_tree["alignment"]["w"] >> w;
    if (realign_tree["alignment"].has_child("zdrop"))
        realign_tree["alignment"]["zdrop"] >> zdrop;
    if (realign_tree["alignment"].has_child("end_bonus"))
        realign_tree["alignment"]["end_bonus"] >> end_bonus;
}

void AlnOpts::print_parameters() {
    std::cerr << "[I::aln::AlnOpts::print_parameters] engine = " << engine << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] flag = " << flag << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] nm_threshold = " << nm_threshold << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] a = " << a << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] b = " << b << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] q = " << q << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] e = " << e << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] q2 = " << q2 << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] e2 = " << e2 << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] w = " << w << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] zdrop = " << zdrop << "\n";
    std::cerr << "[I::aln::AlnOpts::print_parameters] end_bonus = " << end_bonus << "\n";
}

int align_ksw2(
    const char *tseq, const char *qseq, const AlnOpts& opt,
    std::vector<uint32_t>& new_cigar, int& new_score
) {
    int i;
    int8_t a = opt.a, _b = opt.b;
    int8_t b = _b < 0? _b : -_b; // a>0 and b<0
    int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
    int tl = strlen(tseq), ql = strlen(qseq);
    uint8_t *ts, *qs, c[256];
    ksw_extz_t ez;

    memset(&ez, 0, sizeof(ksw_extz_t));
    memset(c, 4, 256);
    c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
    c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
    ts = (uint8_t*)malloc(tl);
    qs = (uint8_t*)malloc(ql);
    for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
    for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
    if (opt.engine == "ksw_extd2_sse") {
        // void ksw_extd2_sse(
        //      void *km, int qlen, const uint8_t *query, int tlen,
        //      const uint8_t *target, int8_t m, const int8_t *mat,
        //      int8_t q, int8_t e, int8_t q2, int8_t e2, int w,
        //      int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
        ksw_extd2_sse(
            0, ql, qs, tl, ts, 5, mat, opt.q, opt.e,
            opt.q2, opt.e2, opt.w, opt.zdrop, opt.end_bonus,
            opt.flag, &ez);
    } else if (opt.engine == "ksw_extz2_sse") {
        // void ksw_extz2_sse(
        //     void *km, int qlen, const uint8_t *query, int tlen, 
        //     const uint8_t *target, int8_t m, const int8_t *mat,
        //     int8_t q, int8_t e, int w, int zdrop,
        //     int end_bonus, int flag, ksw_extz_t *ez)
        ksw_extz2_sse(
            0, ql, qs, tl, ts, 5, mat, opt.q, opt.e, opt.w,
            opt.zdrop, opt.end_bonus, opt.flag, &ez);
    } else {
        // void ksw_extz(
        //     void *km, int qlen, const uint8_t *query, int tlen,
        //     const uint8_t *target, int8_t m, const int8_t *mat,
        //     int8_t gapo, int8_t gape, int w, int zdrop, int flag, ksw_extz_t *ez)
        ksw_extz(
            0, ql, qs, tl, ts, 5, mat, opt.q, opt.e, opt.w,
            opt.zdrop, opt.flag, &ez);
    }
    for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
        // printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
        chain::push_cigar(new_cigar, bam_cigar_oplen(ez.cigar[i]), bam_cigar_op(ez.cigar[i]), false);
    }
    // putchar('\n');
    new_score = ez.score;
    free(ez.cigar); free(ts); free(qs);
    return ez.n_cigar;
}


};
