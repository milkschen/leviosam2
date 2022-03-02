#include "aln.hpp"


namespace Aln {

void AlnOpts::deserialize_realn(ryml::Tree realign_tree) {
    realign_tree["alignment"]["engine"] >> engine;
    realign_tree["alignment"]["flag"] >> flag;
    realign_tree["alignment"]["a"] >> a;
    realign_tree["alignment"]["b"] >> b;
    realign_tree["alignment"]["q"] >> q;
    realign_tree["alignment"]["e"] >> e;
    realign_tree["alignment"]["q2"] >> q2;
    realign_tree["alignment"]["e2"] >> e2;
    realign_tree["alignment"]["w"] >> w;
    realign_tree["alignment"]["zdrop"] >> zdrop;
    realign_tree["alignment"]["end_bonus"] >> end_bonus;
}

void AlnOpts::print_parameters() {
    std::cerr << "engine = " << engine << "\n";
    std::cerr << "flag = " << flag << "\n";
    std::cerr << "a = " << a << "\n";
    std::cerr << "b = " << b << "\n";
    std::cerr << "q = " << q << "\n";
    std::cerr << "e = " << e << "\n";
    std::cerr << "q2 = " << q2 << "\n";
    std::cerr << "e2 = " << e2 << "\n";
    std::cerr << "w = " << w << "\n";
    std::cerr << "zdrop = " << zdrop << "\n";
    std::cerr << "end_bonus = " << end_bonus << "\n";
}

std::vector<uint32_t> align_ksw2(
    const char *tseq, const char *qseq, const AlnOpts& opt
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
    std::vector<uint32_t> new_cigar;
    if (ez.n_cigar == 0)
        std::cerr << "0\n";
    else {
        for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
            printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
            chain::push_cigar(new_cigar, bam_cigar_oplen(ez.cigar[i]), bam_cigar_op(ez.cigar[i]), false);
        }
        putchar('\n');
    }
    free(ez.cigar); free(ts); free(qs);
    return new_cigar;
}


// std::vector<uint32_t> align_extd2(
//     const char *tseq, const char *qseq, const AlnOpts& opt
// ) {
//     int i;
//     int8_t a = opt.a, _b = opt.b;
//     int8_t b = _b < 0? _b : -_b; // a>0 and b<0
//     int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
//     int tl = strlen(tseq), ql = strlen(qseq);
//     uint8_t *ts, *qs, c[256];
//     ksw_extz_t ez;
// 
//     memset(&ez, 0, sizeof(ksw_extz_t));
//     memset(c, 4, 256);
//     c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
//     c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
//     ts = (uint8_t*)malloc(tl);
//     qs = (uint8_t*)malloc(ql);
//     for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
//     for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
//     //ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo, gape, -1, -1, 0, 0, &ez);
//     ksw_extd2_sse(
//         0, ql, qs, tl, ts, 5, mat, opt.q, opt.e,
//         opt.q2, opt.e2, opt.w, opt.zdrop, opt.end_bonus,
//         opt.flag, &ez);
//     // void ksw_extd2_sse(
//     //      void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
//     //      int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
//     std::vector<uint32_t> new_cigar;
//     for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
//         printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//         chain::push_cigar(new_cigar, bam_cigar_oplen(ez.cigar[i]), bam_cigar_op(ez.cigar[i]), false);
//     }
//     putchar('\n');
//     free(ez.cigar); free(ts); free(qs);
//     return new_cigar;
// }

// std::vector<uint32_t> align_extd2(
//     const char *tseq, const char *qseq, int sc_mch, int sc_mis,
//     const int q, const int e, const int8_t q2, const int8_t e2, const int w,
//     const int zdrop, const int end_bonus, const int flag
// )
// {
//     int i;
//     int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
//     int8_t mat[25] = { a,b,b,b,0, b,a,b,b,0, b,b,a,b,0, b,b,b,a,0, 0,0,0,0,0 };
//     int tl = strlen(tseq), ql = strlen(qseq);
//     uint8_t *ts, *qs, c[256];
//     ksw_extz_t ez;
// 
//     memset(&ez, 0, sizeof(ksw_extz_t));
//     memset(c, 4, 256);
//     c['A'] = c['a'] = 0; c['C'] = c['c'] = 1;
//     c['G'] = c['g'] = 2; c['T'] = c['t'] = 3; // build the encoding table
//     ts = (uint8_t*)malloc(tl);
//     qs = (uint8_t*)malloc(ql);
//     for (i = 0; i < tl; ++i) ts[i] = c[(uint8_t)tseq[i]]; // encode to 0/1/2/3
//     for (i = 0; i < ql; ++i) qs[i] = c[(uint8_t)qseq[i]];
//     // ksw_extz(     0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
//     // void ksw_extz(
//     //      void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
//     //      int8_t gapo, int8_t gape, int w, int zdrop, int flag, ksw_extz_t *ez)
//     //ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo, gape, -1, -1, 0, 0, &ez);
//     ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, &ez);
//     // void ksw_extd2_sse(
//     //      void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
//     //      int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
//     std::vector<uint32_t> new_cigar;
//     for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
//         printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
//         chain::push_cigar(new_cigar, bam_cigar_oplen(ez.cigar[i]), bam_cigar_op(ez.cigar[i]), false);
//     }
//     putchar('\n');
//     free(ez.cigar); free(ts); free(qs);
//     return new_cigar;
// }

};
