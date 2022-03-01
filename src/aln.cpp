#include "aln.hpp"

namespace Aln {

std::vector<uint32_t> align_extd2(
    const char *tseq, const char *qseq, int sc_mch, int sc_mis,
    const int q, const int e, const int8_t q2, const int8_t e2, const int w,
    const int zdrop, const int end_bonus, const int flag
)
{
    int i;
    int8_t a = sc_mch, b = sc_mis < 0? sc_mis : -sc_mis; // a>0 and b<0
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
    // ksw_extz(     0, ql, qs, tl, ts, 5, mat, gapo, gape, -1, -1, 0, &ez);
    // void ksw_extz(
    //      void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
    //      int8_t gapo, int8_t gape, int w, int zdrop, int flag, ksw_extz_t *ez)
    //ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, gapo, gape, gapo, gape, -1, -1, 0, 0, &ez);
    ksw_extd2_sse(0, ql, qs, tl, ts, 5, mat, q, e, q2, e2, w, zdrop, end_bonus, flag, &ez);
    // void ksw_extd2_sse(
    //      void *km, int qlen, const uint8_t *query, int tlen, const uint8_t *target, int8_t m, const int8_t *mat,
    //      int8_t q, int8_t e, int8_t q2, int8_t e2, int w, int zdrop, int end_bonus, int flag, ksw_extz_t *ez)
    std::vector<uint32_t> new_cigar;
    for (i = 0; i < ez.n_cigar; ++i) { // print CIGAR
        printf("%d%c", ez.cigar[i]>>4, "MID"[ez.cigar[i]&0xf]);
        chain::push_cigar(new_cigar, bam_cigar_oplen(ez.cigar[i]), bam_cigar_op(ez.cigar[i]), false);
    }
    putchar('\n');
    free(ez.cigar); free(ts); free(qs);
    return new_cigar;
}

};
