/*
 * paffy unanchor: remove the hub (minigraph) anchors of every genome at tandem-repeat loci where the genomes'
 * anchors disagree, so that cactus's BAR realigns those loci from scratch.
 *
 *  Released under the MIT license, see LICENSE.txt
 *
 * Input is one chromosome's PAF as cactus_consolidated reads it (after filter_paf): every target is a node of the
 * hub genome (id=HUB|sN). The minigraph rGFA gives each node of the reference genome its place on the reference
 * path (SN/SO/SR tags), and the node FASTA its bases.
 *
 * Overview:
 * (1) Lay the reference genome's rank-0 nodes out on its path, one contig after another, and confirm each
 *     position with the reference's own PAF records (positions they do not place where the GFA does are
 *     "unframed").
 * (2) Hamming tandem scan of the reference path (periods minPeriod..maxPeriod). Calls merged into long arrays are
 *     satellites and excluded; the rest, padded and merged, are the candidate loci.
 * (3) Walk each non-reference contig's gapless runs in query order. A departure from the reference path (a gap
 *     between two consecutive runs) with both ends in [s-1, e] is inside locus [s, e); one that crosses a locus
 *     boundary, flips strand, goes back or overlaps on the query vetoes it. Each pass of a contig through a locus
 *     is an allele with its own inserted and deleted bases.
 * (4) A locus is cleared when at least minHaps contigs have a pass with both inserted and deleted bases (comp
 *     gate) and nothing vetoes it (span, unframed, crossing, flip, back, qoverlap, shared-alt, string, the SV veto
 *     on bubble loci whose alt nodes hold a non-repeat SV allele, and loci with no alt node and no node boundary).
 * (5) Remove every matched column on a cleared hub interval from every record (the cut), and write the query bases
 *     that lost every hub anchor.
 *
 * The output is byte-identical to the Python reference of the anchor-trim pilots (find_loci.py + cut_paf.py, GFA
 * frame). Every order is a sort or the input order, so the output does not depend on the thread count.
 */

#define _GNU_SOURCE // getline, fseeko, clock_gettime under -std=c99
#include "paf.h"
#include <errno.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <stdarg.h>
#include <zlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Utilities
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void *xrealloc(void *p, size_t n) {
    void *q = realloc(p, n ? n : 1);
    if (q == NULL) {
        st_errAbort("paffy unanchor: out of memory");
    }
    return q;
}

#define VEC(T) struct { T *a; int64_t n, cap; }
#define VEC_GROW(v, k) do { if ((v).n + (k) > (v).cap) { (v).cap = ((v).n + (k)) * 2 + 16; \
        (v).a = xrealloc((v).a, sizeof(*(v).a) * (size_t)(v).cap); } } while (0)
#define VEC_PUSH(v, x) do { VEC_GROW(v, 1); (v).a[(v).n++] = (x); } while (0)
#define VEC_FREE(v) do { free((v).a); (v).a = NULL; (v).n = (v).cap = 0; } while (0)

typedef VEC(int64_t) I64Vec;
typedef VEC(char) CharVec;

static void die2(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    fprintf(stderr, "paffy unanchor ERROR: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(2);
}

static char *vstrprintf(const char *fmt, va_list ap) {
    va_list ap2;
    va_copy(ap2, ap);
    int n = vsnprintf(NULL, 0, fmt, ap2);
    va_end(ap2);
    char *s = st_malloc((size_t)n + 1);
    vsnprintf(s, (size_t)n + 1, fmt, ap);
    return s;
}

static char *strprintf(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *s = vstrprintf(fmt, ap);
    va_end(ap);
    return s;
}

static void cv_printf(CharVec *v, const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *s = vstrprintf(fmt, ap);
    va_end(ap);
    int64_t k = (int64_t)strlen(s);
    VEC_GROW(*v, k + 1);
    memcpy(v->a + v->n, s, (size_t)k + 1);
    v->n += k;
    free(s);
}

static const char *cv_str(CharVec *v) {
    if (v->a == NULL) {
        VEC_GROW(*v, 1);
        v->a[0] = '\0';
    }
    return v->a;
}

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/*
 * The summary log: kept in memory for --logFile and echoed on stderr at log level info.
 */
static stList *log_lines = NULL;

static void ulog(const char *fmt, ...) {
    va_list ap;
    va_start(ap, fmt);
    char *s = vstrprintf(fmt, ap);
    va_end(ap);
    if (log_lines == NULL) {
        log_lines = stList_construct3(0, free);
    }
    stList_append(log_lines, s);
    st_logInfo("%s\n", s);
}

static FILE *open_out(const char *path) {
    FILE *f = fopen(path, "w");
    if (f == NULL) {
        st_errnoAbort("Could not open %s for writing", path);
    }
    return f;
}

static void close_out(FILE *f, const char *path) {
    st_fclose(f, path); // checks the error indicator, then fclose: a short file is never silent
}

static int cmp_i64(const void *a, const void *b) {
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return x < y ? -1 : (x > y ? 1 : 0);
}

// index of the first element > x in sorted a[0..n)
static int64_t upper_bound(const int64_t *a, int64_t n, int64_t x) {
    int64_t lo = 0, hi = n;
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (a[mid] <= x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

// index of the first element >= x in sorted a[0..n)
static int64_t lower_bound(const int64_t *a, int64_t n, int64_t x) {
    int64_t lo = 0, hi = n;
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (a[mid] < x) lo = mid + 1; else hi = mid;
    }
    return lo;
}

static bool starts_with(const char *s, const char *prefix) {
    return strncmp(s, prefix, strlen(prefix)) == 0;
}

static bool is_ascii_space(char c) {
    return c == ' ' || c == '\t' || c == '\n' || c == '\r' || c == '\v' || c == '\f';
}

static bool is_blank_line(const char *s, int64_t n) {
    for (int64_t i = 0; i < n; i++) {
        if (!is_ascii_space(s[i])) return 0;
    }
    return 1;
}

static int64_t parse_int(const char *s, const char *what) {
    char *end;
    errno = 0;
    long long v = strtoll(s, &end, 10);
    if (end == s || *end != '\0' || errno != 0) {
        die2("could not parse %s as an integer: '%s'", what, s);
    }
    return (int64_t)v;
}

// read one line (any length) from a gzFile; returns its length without the newline, or -1 at the end
static int64_t gz_getline(gzFile f, CharVec *buf) {
    buf->n = 0;
    VEC_GROW(*buf, 65536);
    while (1) {
        if (gzgets(f, buf->a + buf->n, (int)(buf->cap - buf->n > INT32_MAX ? INT32_MAX : buf->cap - buf->n)) == NULL) {
            int err;
            const char *msg = gzerror(f, &err);
            if (err != Z_OK && err != Z_STREAM_END) {
                st_errAbort("Error reading a gzipped file: %s", msg);
            }
            if (buf->n == 0) return -1;
            break;
        }
        buf->n += (int64_t)strlen(buf->a + buf->n);
        if (buf->n > 0 && buf->a[buf->n - 1] == '\n') {
            buf->n--;
            buf->a[buf->n] = '\0';
            break;
        }
        VEC_GROW(*buf, buf->cap); // line longer than the buffer: double it and continue
    }
    return buf->n;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Tandem scan (port of the anchor-trim design's hamtr2.c, itself rgfa-collapse -T tandem_cover)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct { int64_t a, b, p; } TRun;
typedef VEC(TRun) TRunVec;

static void scan_period(const char *s, int64_t n, int64_t p, int64_t W, double F, TRunVec *out) {
    int64_t m = n - p, w = p > W ? p : W;
    if (w > m || p >= n) {
        return; // no window fits
    }
    int64_t need = (int64_t)ceil(F * (double)w);
    int64_t hits = 0;
    for (int64_t j = 0; j < w; j++) {
        hits += (s[j] == s[j + p] && s[j] != 'N');
    }
    int64_t cur_s = -1, cur_e = -1;
    for (int64_t i = 0;; i++) { // window starts i in [0, m - w]
        if (hits >= need) {
            int64_t b = i + w + p < n ? i + w + p : n;
            if (cur_s >= 0 && i <= cur_e) {
                if (b > cur_e) cur_e = b;
            } else {
                if (cur_s >= 0) {
                    TRun r = { cur_s, cur_e, p };
                    VEC_PUSH(*out, r);
                }
                cur_s = i;
                cur_e = b;
            }
        }
        if (i + w >= m) {
            break;
        }
        hits += (s[i + w] == s[i + w + p] && s[i + w] != 'N') - (s[i] == s[i + p] && s[i] != 'N');
    }
    if (cur_s >= 0) {
        TRun r = { cur_s, cur_e, p };
        VEC_PUSH(*out, r);
    }
}

static int cmp_trun_ab(const void *x, const void *y) {
    const TRun *a = x, *b = y;
    if (a->a != b->a) return a->a < b->a ? -1 : 1;
    if (a->b != b->b) return a->b < b->b ? -1 : 1;
    return a->p < b->p ? -1 : (a->p > b->p ? 1 : 0);
}

int64_t *tandem_scan(const char *seq, int64_t n, int64_t min_period, int64_t max_period, int64_t window,
                     double min_identity, int64_t threads, int64_t *call_number) {
    if (max_period > 100) {
        st_errAbort("tandem_scan: the maximum period must be at most 100");
    }
    int64_t np = max_period >= min_period ? max_period - min_period + 1 : 0;
    TRunVec *per = st_calloc(np > 0 ? np : 1, sizeof(TRunVec));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(threads > 0 ? threads : 1)
#endif
    for (int64_t k = 0; k < np; k++) {
        scan_period(seq, n, min_period + k, window, min_identity, &per[k]);
    }
    (void)threads;
    TRunVec all = { 0 };
    for (int64_t k = 0; k < np; k++) {
        if (per[k].n > 0) {
            VEC_GROW(all, per[k].n);
            memcpy(all.a + all.n, per[k].a, sizeof(TRun) * (size_t)per[k].n);
            all.n += per[k].n;
        }
        VEC_FREE(per[k]);
    }
    free(per);
    if (all.n > 0) {
        qsort(all.a, (size_t)all.n, sizeof(TRun), cmp_trun_ab);
    }
    // union runs (touching runs are one), each labelled with the smallest period whose covered bp, summed over
    // that period's runs starting in the union run, is within 10% of the best period's
    I64Vec calls = { 0 };
    int64_t tally[101];
    int64_t i = 0;
    while (i < all.n) {
        int64_t us = all.a[i].a, ue = all.a[i].b, j = i + 1;
        while (j < all.n && all.a[j].a <= ue) {
            if (all.a[j].b > ue) ue = all.a[j].b;
            j++;
        }
        memset(tally, 0, sizeof(tally));
        for (int64_t k = i; k < j; k++) {
            tally[all.a[k].p] += all.a[k].b - all.a[k].a;
        }
        int64_t bb = -1;
        for (int64_t k = i; k < j; k++) {
            if (tally[all.a[k].p] > bb) bb = tally[all.a[k].p];
        }
        int64_t lab = 0;
        for (int64_t p = min_period; p <= max_period && p <= 100; p++) {
            if (tally[p] * 10 >= bb * 9) {
                lab = p;
                break;
            }
        }
        VEC_PUSH(calls, us);
        VEC_PUSH(calls, ue);
        VEC_PUSH(calls, lab);
        i = j;
    }
    VEC_FREE(all);
    *call_number = calls.n / 3;
    return calls.a;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Parameters
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum gate_kind { GATE_COMP, GATE_VARIABLE, GATE_NONE };

typedef struct {
    char *paf, *nodes, *gfa, *ref, *hub;
    char *out_paf, *loci_bed, *hub_bed, *veto_bed, *query_bed, *calls_bed, *sv_report, *log_file, *scan_fasta;
    bool dry_run, sv_veto, keep_cigar_no_boundary;
    int gate;
    int64_t min_haps, max_len, max_span, max_string, caf_trim, min_period, max_period, window, pad, pad_periods;
    int64_t merge, array_merge, array_min, sv_min_bp, sv_k, sv_tight, sv_wide, threads;
    double min_ident, sv_max_local, sv_inv_min, sv_inv_margin;
} Params;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Hub nodes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum {
    NODE_RANK0 = 1, NODE_GFA_ALT = 2, NODE_FASTA = 4, NODE_ALT_OUT = 8, NODE_REF_ALT = 16, NODE_SHARED = 32,
    NODE_CLEARED_ALT = 64
};

typedef struct {
    char *name;           // with the id=HUB| prefix
    char *rank0_contig;   // SN contig of a rank-0 node
    int64_t so, ln;       // GFA SO and LN (ln -1 if not in the GFA)
    int64_t base;         // global position of node offset 0 (rank-0 nodes in the layout), else -1
    int64_t gfa_order;    // order of the node's first S line
    char *seq;            // node FASTA bases (uppercase), NULL if missing
    int64_t seq_len;
    int64_t tlen;         // PAF column 7 of the last record targeting it, -1 if none
    int64_t rank;         // byte order of the name among all nodes
    int64_t clear_locus;  // locus that clears the whole (alt) node, -1 if none
    int flags;
} Node;

typedef struct {
    stHash *index;        // name -> index + 1
    VEC(Node) v;
} NodeTable;

static int64_t node_get(NodeTable *nt, const char *name, bool create) {
    void *x = stHash_search(nt->index, (void *)name);
    if (x != NULL) {
        return (int64_t)(intptr_t)x - 1;
    }
    if (!create) {
        return -1;
    }
    Node nd;
    memset(&nd, 0, sizeof(nd));
    nd.name = stString_copy(name);
    nd.ln = -1;
    nd.base = -1;
    nd.so = -1;
    nd.tlen = -1;
    nd.gfa_order = -1;
    nd.clear_locus = -1;
    VEC_PUSH(nt->v, nd);
    stHash_insert(nt->index, nd.name, (void *)(intptr_t)nt->v.n);
    return nt->v.n - 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// PAF records
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct {
    char **f;
    int64_t n, cap;
} Fields;

static void split_fields(char *line, Fields *fl) {
    fl->n = 0;
    char *p = line;
    while (1) {
        if (fl->n == fl->cap) {
            fl->cap = fl->cap * 2 + 32;
            fl->f = xrealloc(fl->f, sizeof(char *) * (size_t)fl->cap);
        }
        fl->f[fl->n++] = p;
        char *t = strchr(p, '\t');
        if (t == NULL) break;
        *t = '\0';
        p = t + 1;
    }
}

typedef struct { int64_t len; char op; } Op;
typedef VEC(Op) OpVec;

typedef struct {
    char *q, *t;
    int64_t qlen, qs, qe, tlen, ts, te, nmatch;
    bool rev, secondary;
    int64_t cg_idx; // field index of the first cg:Z: tag
} Rec;

static void parse_cigar(const char *s, OpVec *ops) {
    ops->n = 0;
    const char *p = s;
    while (*p != '\0') {
        if (*p < '0' || *p > '9') {
            die2("cigar op without a length in cg:Z:%s", s);
        }
        int64_t len = 0;
        while (*p >= '0' && *p <= '9') {
            len = len * 10 + (*p - '0');
            p++;
        }
        char op = *p;
        if (op != 'M' && op != 'I' && op != 'D' && op != 'X' && op != '=') {
            die2("cigar op '%c' is not one of MIDX= in cg:Z:%s", op == '\0' ? ' ' : op, s);
        }
        p++;
        Op o = { len, op };
        VEC_PUSH(*ops, o);
    }
}

// parse a split PAF line (fields modified in place) into r, and its cigar into ops
static void parse_record(Fields *fl, Rec *r, OpVec *ops) {
    if (fl->n < 12) {
        die2("PAF record with %" PRIi64 " fields (at least 12 needed): %s", fl->n, fl->f[0]);
    }
    char **t = fl->f;
    r->q = t[0];
    r->qlen = parse_int(t[1], "the query length");
    r->qs = parse_int(t[2], "the query start");
    r->qe = parse_int(t[3], "the query end");
    if (strcmp(t[4], "+") == 0) {
        r->rev = 0;
    } else if (strcmp(t[4], "-") == 0) {
        r->rev = 1;
    } else {
        die2("PAF record with strand '%s': %s", t[4], t[0]);
    }
    r->t = t[5];
    r->tlen = parse_int(t[6], "the target length");
    r->ts = parse_int(t[7], "the target start");
    r->te = parse_int(t[8], "the target end");
    r->nmatch = parse_int(t[9], "the number of matches");
    r->secondary = 0;
    r->cg_idx = -1;
    for (int64_t k = 12; k < fl->n; k++) {
        if (r->cg_idx < 0 && starts_with(t[k], "cg:Z:")) {
            r->cg_idx = k;
        }
        if (strcmp(t[k], "tp:A:S") == 0) {
            r->secondary = 1;
        }
    }
    if (r->cg_idx < 0) {
        die2("PAF record without a cg:Z: tag: %s %s", t[0], t[5]);
    }
    if (ops != NULL) {
        parse_cigar(t[r->cg_idx] + 5, ops);
    }
}

typedef struct { int64_t qf, n, t0; } Run;
typedef VEC(Run) RunVec;

// gapless runs of a record (neighbouring =/X/M ops merged): qf is the run's lowest forward-strand query base
static void record_runs(Rec *r, OpVec *ops, RunVec *runs) {
    runs->n = 0;
    int64_t q = 0, tt = r->ts, last_q = -1, last_t = -1;
    for (int64_t k = 0; k < ops->n; k++) {
        int64_t n = ops->a[k].len;
        char op = ops->a[k].op;
        if (op == 'M' || op == '=' || op == 'X') {
            if (runs->n > 0 && last_q == q && last_t == tt) {
                runs->a[runs->n - 1].n += n;
            } else {
                Run x = { q, n, tt };
                VEC_PUSH(*runs, x);
            }
            q += n;
            tt += n;
            last_q = q;
            last_t = tt;
        } else if (op == 'I') {
            q += n;
        } else {
            tt += n;
        }
    }
    for (int64_t k = 0; k < runs->n; k++) {
        Run *x = &runs->a[k];
        x->qf = r->rev ? r->qe - x->qf - x->n : r->qs + x->qf;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The state of one run
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum { R_SPAN, R_UNFRAMED, R_CROSSING, R_FLIP, R_BACK, R_QOVERLAP, R_SHARED_ALT, R_STRING, R_GATE,
       R_SV_ORIENT, R_SV_INV, R_SV_NOVEL, R_CIGAR_NOBOUNDARY, R_NUM };
static const char *reason_names[R_NUM] = { "span", "unframed", "crossing", "flip", "back", "qoverlap", "shared-alt",
                                           "string", "gate", "sv-orient", "sv-inv", "sv-novel", "cigar-noboundary" };

enum { K_BACK, K_FLIP, K_GAP, K_QOVERLAP, K_NUM };
static const char *kind_names[K_NUM] = { "back", "flip", "gap", "qoverlap" };

typedef struct {
    int64_t s, e, p;
} Locus;

typedef struct {
    int64_t li, contig, passno, ins, del, alt0, nalt, order;
    int dir;
    bool has_alt, has_replace;
} Pass;

typedef struct {
    int64_t node, ri;
    bool rev;
} AltUse;

typedef struct { int64_t node, li; } AltIn;

typedef struct {
    char *name;
    int64_t genome;            // index of the contig's genome (for nHaps / genomes)
    VEC(int64_t) recs;         // (file offset, record index) pairs of its primary records
} QContig;

typedef struct {
    int64_t qf, n, rlo, rhi, ri, order, node;
    int o;                     // +1 / -1 for runs on rank-0 nodes, 0 for runs on alt nodes
    bool rev;                  // the record's strand
} DRun;

typedef VEC(DRun) DRunVec;

typedef struct { int64_t li, node, ts, te; bool alt; } HubIv;

typedef struct {
    Params *a;
    NodeTable nt;
    char *ref_q, *hub_t;
    // reference layout
    int64_t ncontigs;
    char **contigs;
    int64_t *clen, *coff;
    int64_t G;
    char *seq;
    I64Vec bounds;
    I64Vec uf_s, uf_emax;     // unframed intervals: starts and running max of ends
    // loci
    VEC(Locus) loci;
    I64Vec lstarts;
    // query contigs (non-reference primaries)
    stHash *qindex;
    VEC(QContig) qc;
    stHash *gindex;
    int64_t ngenomes;
    // departures
    int64_t *inside, *contig_end_inside, *maxstr, *minstr, *comp, *mech; // mech: 3 per locus
    uint32_t *mask;
    bool *bubble;
    VEC(Pass) passes;
    VEC(AltUse) altuses;
    VEC(AltIn) altin;
    int64_t nkind[K_NUM];
    CharVec *sv_detail;       // per locus (allocated lazily)
} State;

static int64_t locus_of(State *st, int64_t x) {
    int64_t i = upper_bound(st->lstarts.a, st->lstarts.n, x) - 1;
    return (i >= 0 && x < st->loci.a[i].e) ? i : -1;
}

static int64_t genome_index(State *st, const char *q) {
    // genome_of: q[3:].split('|')[0] if q starts with id= else q.split('|')[0]
    const char *s = starts_with(q, "id=") ? q + 3 : q;
    const char *bar = strchr(s, '|');
    char *g = bar == NULL ? stString_copy(s) : stString_getSubString(s, 0, bar - s);
    void *x = stHash_search(st->gindex, g);
    if (x != NULL) {
        free(g);
        return (int64_t)(intptr_t)x - 1;
    }
    stHash_insert(st->gindex, g, (void *)(intptr_t)(++st->ngenomes));
    return st->ngenomes - 1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Inputs: GFA S lines and the node FASTA
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static gzFile open_gz(const char *path) {
    gzFile f = gzopen(path, "r");
    if (f == NULL) {
        st_errnoAbort("Could not open %s", path);
    }
    gzbuffer(f, 1 << 20);
    return f;
}

static void read_gfa(State *st, int64_t *n_s, int64_t *n_rank0, int64_t *n_alt) {
    Params *a = st->a;
    gzFile f = open_gz(a->gfa);
    CharVec line = { 0 };
    Fields fl = { 0 };
    char *ref_hash = strprintf("%s#", a->ref);
    // cactus event_to_pansn_prefix: a reference named S.N (N all digits) is S#N# in PanSN
    char *ref_pansn = NULL;
    char *dot = strrchr(a->ref, '.');
    if (dot != NULL && dot[1] != '\0' && strspn(dot + 1, "0123456789") == strlen(dot + 1)) {
        ref_pansn = strprintf("%.*s#%s#", (int)(dot - a->ref), a->ref, dot + 1);
    }
    *n_s = *n_rank0 = *n_alt = 0;
    while (gz_getline(f, &line) >= 0) {
        if (line.n < 2 || line.a[0] != 'S' || line.a[1] != '\t') {
            continue;
        }
        (*n_s)++;
        split_fields(line.a, &fl);
        if (fl.n < 3) {
            die2("GFA S line with fewer than 3 fields");
        }
        char *name = strprintf("%s%s", st->hub_t, fl.f[1]);
        const char *ln_s = NULL, *sn = NULL, *so_s = NULL, *sr_s = NULL;
        for (int64_t k = 3; k < fl.n; k++) { // later tags override earlier ones
            const char *tg = fl.f[k];
            const char *val = strlen(tg) >= 5 ? tg + 5 : "";
            if (strncmp(tg, "LN", 2) == 0) ln_s = val;
            else if (strncmp(tg, "SN", 2) == 0) sn = val;
            else if (strncmp(tg, "SO", 2) == 0) so_s = val;
            else if (strncmp(tg, "SR", 2) == 0) sr_s = val;
        }
        int64_t ln = ln_s != NULL ? parse_int(ln_s, "a GFA LN tag") : (int64_t)strlen(fl.f[2]);
        int64_t so = so_s != NULL ? parse_int(so_s, "a GFA SO tag") : -1;
        int64_t sr = sr_s != NULL ? parse_int(sr_s, "a GFA SR tag") : -1;
        char *contig = NULL;
        if (sn != NULL && sr == 0) {
            if (starts_with(sn, st->ref_q)) {
                contig = stString_copy(sn + strlen(st->ref_q));
            } else if (starts_with(sn, ref_hash) || (ref_pansn != NULL && starts_with(sn, ref_pansn))) {
                const char *h1 = strchr(sn, '#');
                const char *h2 = h1 == NULL ? NULL : strchr(h1 + 1, '#');
                if (h2 == NULL) {
                    die2("GFA SN tag %s has no contig field", sn);
                }
                contig = stString_copy(h2 + 1);
            }
        }
        int64_t i = node_get(&st->nt, name, 1);
        free(name);
        Node *nd = &st->nt.v.a[i];
        if (nd->gfa_order < 0) {
            nd->gfa_order = *n_s - 1;
        }
        nd->ln = ln;
        if (contig != NULL) {
            if (!(nd->flags & NODE_RANK0)) (*n_rank0)++;
            nd->flags |= NODE_RANK0;
            free(nd->rank0_contig);
            nd->rank0_contig = contig;
            nd->so = so;
        } else {
            if (!(nd->flags & NODE_GFA_ALT)) (*n_alt)++;
            nd->flags |= NODE_GFA_ALT;
        }
    }
    gzclose(f);
    VEC_FREE(line);
    free(fl.f);
    free(ref_hash);
    free(ref_pansn);
}

typedef struct { char *name; char *seq; int64_t len; } FaRec;
typedef VEC(FaRec) FaVec;

// read a FASTA (plain or gzip): names are the first word of the header, bases uppercased (ASCII only); a repeated
// name keeps its first place and its last sequence
static void read_fasta(const char *path, FaVec *out) {
    gzFile f = open_gz(path);
    CharVec line = { 0 }, seq = { 0 };
    stHash *seen = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    char *name = NULL;
    while (1) {
        int64_t n = gz_getline(f, &line);
        if (n < 0 || (n > 0 && line.a[0] == '>')) {
            if (name != NULL) {
                for (int64_t i = 0; i < seq.n; i++) {
                    if (seq.a[i] >= 'a' && seq.a[i] <= 'z') seq.a[i] -= 32;
                }
                char *s = st_malloc((size_t)seq.n + 1);
                if (seq.n > 0) memcpy(s, seq.a, (size_t)seq.n); // seq.a is NULL before any sequence line
                s[seq.n] = '\0';
                void *x = stHash_search(seen, name);
                if (x != NULL) {
                    FaRec *r = &out->a[(int64_t)(intptr_t)x - 1];
                    free(r->seq);
                    r->seq = s;
                    r->len = seq.n;
                    free(name);
                } else {
                    FaRec r = { name, s, seq.n };
                    VEC_PUSH(*out, r);
                    stHash_insert(seen, name, (void *)(intptr_t)out->n);
                }
                name = NULL;
            }
            if (n < 0) break;
            const char *h = line.a + 1;
            while (*h != '\0' && is_ascii_space(*h)) h++;
            const char *e = h;
            while (*e != '\0' && !is_ascii_space(*e)) e++;
            if (e == h) {
                die2("FASTA header without a name in %s", path);
            }
            name = stString_getSubString(h, 0, e - h);
            seq.n = 0;
        } else if (name != NULL) {
            int64_t b = 0, e = n;
            while (b < e && is_ascii_space(line.a[b])) b++;
            while (e > b && is_ascii_space(line.a[e - 1])) e--;
            VEC_GROW(seq, e - b);
            if (e > b) memcpy(seq.a + seq.n, line.a + b, (size_t)(e - b));
            seq.n += e - b;
        }
    }
    gzclose(f);
    stHash_destruct(seen);
    VEC_FREE(line);
    VEC_FREE(seq);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Scan only (--scanFasta)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static int scan_only(Params *a) {
    double t0 = now_seconds();
    FaVec fa = { 0 };
    read_fasta(a->scan_fasta, &fa);
    int64_t **calls = st_calloc(fa.n > 0 ? fa.n : 1, sizeof(int64_t *));
    int64_t *ncalls = st_calloc(fa.n > 0 ? fa.n : 1, sizeof(int64_t));
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1) num_threads(a->threads > 0 ? a->threads : 1)
#endif
    for (int64_t i = 0; i < fa.n; i++) {
        calls[i] = tandem_scan(fa.a[i].seq, fa.a[i].len, a->min_period, a->max_period, a->window, a->min_ident, 1,
                               &ncalls[i]);
    }
    FILE *f = open_out(a->calls_bed);
    for (int64_t i = 0; i < fa.n; i++) {
        int64_t bp = 0;
        for (int64_t k = 0; k < ncalls[i]; k++) {
            int64_t *c = calls[i] + 3 * k;
            fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n", fa.a[i].name, c[0], c[1], c[2]);
            bp += c[1] - c[0];
        }
        ulog("%s: %" PRIi64 " calls, %" PRIi64 " bp", fa.a[i].name, ncalls[i], bp);
        free(calls[i]);
        free(fa.a[i].name);
        free(fa.a[i].seq);
    }
    close_out(f, a->calls_bed);
    ulog("scan time %.1fs", now_seconds() - t0);
    free(calls);
    free(ncalls);
    VEC_FREE(fa);
    return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Departures
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static int cmp_drun(const void *x, const void *y) {
    const DRun *a = x, *b = y;
    if (a->qf != b->qf) return a->qf < b->qf ? -1 : 1;
    if (a->n != b->n) return a->n < b->n ? -1 : 1;
    if (a->node != b->node) return a->node < b->node ? -1 : 1; // node holds the name rank here
    int64_t ka = a->o ? a->rlo : -1, kb = b->o ? b->rlo : -1;
    if (ka != kb) return ka < kb ? -1 : 1;
    return a->order < b->order ? -1 : (a->order > b->order ? 1 : 0);
}

typedef struct {
    CharVec line;
    Fields fl;
    OpVec ops;
    RunVec runs;
} Scratch;

static int64_t read_line_at(FILE *f, int64_t offset, CharVec *line) {
    if (fseeko(f, (off_t)offset, SEEK_SET) != 0) {
        st_errnoAbort("Could not seek in the input PAF");
    }
    char *buf = line->a;
    size_t cap = (size_t)line->cap;
    ssize_t n = getline(&buf, &cap, f);
    line->a = buf;
    line->cap = (int64_t)cap;
    if (n < 0) {
        st_errAbort("Could not re-read a record of the input PAF: was it changed while running?");
    }
    if (n > 0 && buf[n - 1] == '\n') buf[--n] = '\0';
    line->n = n;
    return n;
}

// the gapless runs of one query contig's primary records, sorted as the walk takes them
static void load_contig_runs(State *st, FILE *in, QContig *qc, Scratch *sc, DRunVec *out) {
    out->n = 0;
    int64_t order = 0;
    for (int64_t k = 0; k < qc->recs.n; k += 2) {
        int64_t ri = qc->recs.a[k + 1];
        read_line_at(in, qc->recs.a[k], &sc->line);
        split_fields(sc->line.a, &sc->fl);
        Rec r;
        parse_record(&sc->fl, &r, &sc->ops);
        record_runs(&r, &sc->ops, &sc->runs);
        int64_t node = node_get(&st->nt, r.t, 0);
        assert(node >= 0);
        Node *nd = &st->nt.v.a[node];
        for (int64_t j = 0; j < sc->runs.n; j++) {
            Run *x = &sc->runs.a[j];
            DRun d;
            d.qf = x->qf;
            d.n = x->n;
            d.ri = ri;
            d.order = order++;
            d.node = node;
            d.rev = r.rev;
            if (nd->base >= 0) {
                d.rlo = nd->base + x->t0;
                d.rhi = nd->base + x->t0 + x->n - 1;
                d.o = r.rev ? -1 : 1;
            } else {
                d.rlo = d.rhi = -1;
                d.o = 0;
            }
            VEC_PUSH(*out, d);
        }
    }
    // sort by (qf, n, node name, rlo or -1), ties in insertion order
    for (int64_t j = 0; j < out->n; j++) out->a[j].node = st->nt.v.a[out->a[j].node].rank;
    if (out->n > 0) qsort(out->a, (size_t)out->n, sizeof(DRun), cmp_drun);
}

static void mark_alt_out(State *st, AltUse *alts, int64_t n) {
    for (int64_t k = 0; k < n; k++) st->nt.v.a[alts[k].node].flags |= NODE_ALT_OUT;
}

static void walk_contig(State *st, int64_t ci, DRunVec *runs, int64_t *rank_to_node, int64_t *passno_stamp,
                        int64_t *passno_val, bool *multi_pass_locus, int64_t *npass2) {
    Locus *L = st->loci.a;
    int64_t nl = st->loci.n;
    // contig ends inside a locus
    int64_t first_p = -1, last_p = -1;
    for (int64_t k = 0; k < runs->n; k++) {
        if (runs->a[k].o != 0) {
            if (first_p < 0) first_p = k;
            last_p = k;
        }
    }
    if (first_p >= 0) {
        for (int side = 0; side < 2; side++) {
            int64_t k = side == 0 ? first_p : last_p;
            DRun *r = &runs->a[k];
            int64_t x = ((side == 0) == (r->o == 1)) ? r->rlo : r->rhi;
            int64_t li = locus_of(st, x);
            if (li >= 0 && (side == 0 ? k == 0 : k == runs->n - 1)) {
                st->contig_end_inside[li]++;
            }
        }
    }
    VEC(AltUse) alts = { 0 };
    DRun *prev = NULL;
    int64_t last_in = -1, cur_pass = -1;
    for (int64_t k = 0; k < runs->n; k++) {
        DRun *r = &runs->a[k];
        if (r->o == 0) {
            AltUse u = { rank_to_node[r->node], r->ri, r->rev };
            VEC_PUSH(alts, u);
            continue;
        }
        if (prev == NULL) {
            mark_alt_out(st, alts.a, alts.n);
        } else {
            int64_t dq = r->qf - (prev->qf + prev->n), dr = 0, x, y;
            int kind;
            if (prev->o != r->o) {
                kind = K_FLIP;
                if (prev->o == 1) { x = prev->rhi; y = r->rlo; } else { x = prev->rlo; y = r->rhi; }
                if (x > y) { int64_t t = x; x = y; y = t; }
            } else if (r->o == 1) {
                x = prev->rhi; y = r->rlo; dr = r->rlo - prev->rhi - 1; kind = K_GAP;
            } else {
                x = r->rhi; y = prev->rlo; dr = prev->rlo - r->rhi - 1; kind = K_GAP;
            }
            if (kind == K_GAP && dr < 0) {
                kind = K_BACK;
                if (x > y) { int64_t t = x; x = y; y = t; }
            }
            if (kind == K_GAP && dq < 0) {
                kind = K_QOVERLAP;
            }
            if (!(kind == K_GAP && dq == 0 && dr == 0)) {
                st->nkind[kind]++;
                int64_t last_in_prev = last_in;
                last_in = -1;
                if (kind == K_GAP) {
                    int64_t ins = -1;
                    // loci with s <= y and e > x, in descending order
                    for (int64_t li = upper_bound(st->lstarts.a, nl, y) - 1; li >= 0 && L[li].e > x; li--) {
                        if (L[li].s - 1 <= x && y <= L[li].e) {
                            ins = li;
                        } else if (x < L[li].s - 1 && y > L[li].e) {
                            // spans the locus: its adjacency passes over it
                        } else {
                            st->mask[li] |= 1u << R_CROSSING;
                        }
                    }
                    if (ins >= 0) {
                        st->inside[ins]++;
                        if (last_in_prev != ins) {
                            if (passno_stamp[ins] != ci) {
                                passno_stamp[ins] = ci;
                                passno_val[ins] = 0;
                            }
                            passno_val[ins]++;
                            if (passno_val[ins] == 2) {
                                (*npass2)++;
                                multi_pass_locus[ins] = 1;
                            }
                            Pass p;
                            memset(&p, 0, sizeof(p));
                            p.li = ins;
                            p.contig = ci;
                            p.passno = passno_val[ins];
                            p.dir = r->o;
                            p.alt0 = st->altuses.n;
                            p.order = st->passes.n;
                            VEC_PUSH(st->passes, p);
                            cur_pass = st->passes.n - 1;
                        }
                        last_in = ins;
                        Pass *p = &st->passes.a[cur_pass];
                        assert(p->li == ins && p->alt0 + p->nalt == st->altuses.n);
                        p->ins += dq;
                        p->del += dr;
                        if (alts.n > 0) p->has_alt = 1;
                        if (dq > 0 && dr > 0) p->has_replace = 1;
                        for (int64_t j = 0; j < alts.n; j++) {
                            VEC_PUSH(st->altuses, alts.a[j]);
                            AltIn ai = { alts.a[j].node, ins };
                            VEC_PUSH(st->altin, ai);
                        }
                        p->nalt += alts.n;
                    } else {
                        mark_alt_out(st, alts.a, alts.n);
                    }
                } else {
                    int64_t lx = locus_of(st, x), ly = locus_of(st, y);
                    if (lx >= 0 && !(x < L[lx].s - 1 && y > L[lx].e)) st->mask[lx] |= 1u << (kind == K_FLIP ? R_FLIP : kind == K_BACK ? R_BACK : R_QOVERLAP);
                    if (ly >= 0 && ly != lx && !(x < L[ly].s - 1 && y > L[ly].e)) st->mask[ly] |= 1u << (kind == K_FLIP ? R_FLIP : kind == K_BACK ? R_BACK : R_QOVERLAP);
                    mark_alt_out(st, alts.a, alts.n);
                }
            }
            // on a trivial step the alt runs since prev are discarded
        }
        prev = r;
        alts.n = 0;
    }
    mark_alt_out(st, alts.a, alts.n); // alt runs after the last anchored run
    VEC_FREE(alts);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// SV features (k-mer containment; k <= 16)
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct { uint64_t hi, lo; } Kmer;
typedef VEC(Kmer) KmerVec;

static int cmp_kmer(const void *x, const void *y) {
    const Kmer *a = x, *b = y;
    if (a->hi != b->hi) return a->hi < b->hi ? -1 : 1;
    return a->lo < b->lo ? -1 : (a->lo > b->lo ? 1 : 0);
}

// the distinct k-mers of s[0..n) with no N, sorted
static void kset(const char *s, int64_t n, int64_t K, KmerVec *out) {
    out->n = 0;
    int64_t last_n = -1; // last N position
    for (int64_t i = 0; i < n; i++) {
        if (s[i] == 'N') last_n = i;
        int64_t st0 = i - K + 1;
        if (st0 >= 0 && last_n < st0) {
            unsigned char buf[16] = { 0 };
            memcpy(buf, s + st0, (size_t)K);
            Kmer km;
            memcpy(&km.hi, buf, 8);
            memcpy(&km.lo, buf + 8, 8);
            VEC_PUSH(*out, km);
        }
    }
    if (out->n > 0) {
        qsort(out->a, (size_t)out->n, sizeof(Kmer), cmp_kmer);
        int64_t j = 0;
        for (int64_t i = 0; i < out->n; i++) {
            if (j == 0 || cmp_kmer(&out->a[j - 1], &out->a[i]) != 0) out->a[j++] = out->a[i];
        }
        out->n = j;
    }
}

static double contain(KmerVec *q, KmerVec *t) {
    if (q->n == 0) return 0.0;
    int64_t i = 0, j = 0, hit = 0;
    while (i < q->n && j < t->n) {
        int c = cmp_kmer(&q->a[i], &t->a[j]);
        if (c == 0) { hit++; i++; j++; }
        else if (c < 0) i++;
        else j++;
    }
    return (double)hit / (double)q->n;
}

static char *revcomp(const char *s, int64_t n) {
    char *r = st_malloc((size_t)n + 1);
    for (int64_t i = 0; i < n; i++) {
        char c = s[n - 1 - i];
        switch (c) {
            case 'A': c = 'T'; break;
            case 'C': c = 'G'; break;
            case 'G': c = 'C'; break;
            case 'T': c = 'A'; break;
            default: break; // N and every other byte unchanged
        }
        r[i] = c;
    }
    r[n] = '\0';
    return r;
}

typedef struct { int64_t rank, node, pass, ri, contig; bool fwd; } SvUse;

static int cmp_svuse(const void *x, const void *y) {
    const SvUse *a = x, *b = y;
    if (a->rank != b->rank) return a->rank < b->rank ? -1 : 1;
    if (a->pass != b->pass) return a->pass < b->pass ? -1 : 1;
    return a->ri < b->ri ? -1 : (a->ri > b->ri ? 1 : 0);
}

typedef struct {
    int64_t li, len, fwd, rev, ncontigs, ngenomes, nontr;
    char *name;
    double kf, kfw, krw;
    uint32_t why;
} SvRow;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The cut
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef struct { int64_t *s, *e, n; } Cleared;

typedef struct { int64_t q0, t0; OpVec ops; } Frag;
typedef struct { Frag *a; int64_t n, cap; } FragVec; // entries up to cap are zeroed or own their ops

static Frag *frag_open(FragVec *v) {
    if (v->n == v->cap) {
        int64_t cap = v->cap * 2 + 8;
        v->a = xrealloc(v->a, sizeof(Frag) * (size_t)cap);
        memset(v->a + v->cap, 0, sizeof(Frag) * (size_t)(cap - v->cap));
        v->cap = cap;
    }
    Frag *f = &v->a[v->n];
    f->ops.n = 0;
    return f;
}

static void push_op(OpVec *v, int64_t len, char op) {
    if (v->n > 0 && v->a[v->n - 1].op == op) {
        v->a[v->n - 1].len += len;
    } else {
        Op o = { len, op };
        VEC_PUSH(*v, o);
    }
}

typedef struct { int64_t s, e; } Iv;
typedef VEC(Iv) IvVec;

static int cmp_iv(const void *x, const void *y) {
    const Iv *a = x, *b = y;
    if (a->s != b->s) return a->s < b->s ? -1 : 1;
    return a->e < b->e ? -1 : (a->e > b->e ? 1 : 0);
}

static void merge_ivs(IvVec *v) {
    if (v->n == 0) return;
    qsort(v->a, (size_t)v->n, sizeof(Iv), cmp_iv);
    int64_t j = 0;
    for (int64_t i = 0; i < v->n; i++) {
        if (j > 0 && v->a[i].s <= v->a[j - 1].e) {
            if (v->a[i].e > v->a[j - 1].e) v->a[j - 1].e = v->a[i].e;
        } else {
            v->a[j++] = v->a[i];
        }
    }
    v->n = j;
}

// a satellite array A with A.start <= e and A.end > s (the last array starting at or before e)
static bool in_array(I64Vec *astarts, IvVec *arrays, int64_t s, int64_t e) {
    int64_t i = upper_bound(astarts->a, astarts->n, e) - 1;
    return i >= 0 && arrays->a[i].e > s;
}

// add a record's matched query interval to v, extending the last one when they touch: a cigar's =/X ops (and the
// pieces the cut makes of them) are contiguous on the query, upward on + and downward on -, so this keeps v to about
// one interval per gapless stretch instead of one per op. Only the union of v is ever used (merge_ivs at the end).
static void add_query_iv(IvVec *v, Rec *r, int64_t q, int64_t n) {
    Iv x;
    if (r->rev) { x.s = r->qe - q - n; x.e = r->qe - q; } else { x.s = r->qs + q; x.e = r->qs + q + n; }
    if (v->n > 0) {
        Iv *l = &v->a[v->n - 1];
        if (x.s == l->e) { l->e = x.e; return; }
        if (x.e == l->s) { l->s = x.s; return; }
    }
    VEC_PUSH(*v, x);
}

// cut_paf.py cut_record: returns 0 if no matched column is on a cleared position (the record is unchanged), else
// 1 with the fragments in frags (possibly none). Adds the record's matched query intervals to removed and the
// surviving ones to kept.
static bool cut_record(Rec *r, OpVec *ops, Cleared *c, FragVec *frags, IvVec *removed, IvVec *kept) {
    int64_t j0 = upper_bound(c->e, c->n, r->ts);
    if (j0 >= c->n || c->s[j0] >= r->te) {
        return 0;
    }
    int64_t q = 0, tt = r->ts;
    frags->n = 0;
    bool removed_any = 0, open = 0;
    OpVec pending = { 0 };
    int64_t kept0 = kept->n; // the caller passes kept empty: add_query_iv may extend its last interval
    for (int64_t k = 0; k < ops->n; k++) {
        int64_t n = ops->a[k].len;
        char op = ops->a[k].op;
        if (op == 'M' || op == '=' || op == 'X') {
            int64_t pos = tt, end = tt + n;
            int64_t j = upper_bound(c->e, c->n, pos); // first interval ending after pos
            while (pos < end) {
                int64_t stop;
                bool clr;
                if (j < c->n && c->s[j] <= pos && pos < c->e[j]) {
                    stop = end < c->e[j] ? end : c->e[j];
                    clr = 1;
                    j++;
                } else {
                    int64_t nxt = j < c->n ? c->s[j] : end;
                    stop = end < nxt ? end : nxt;
                    clr = 0;
                }
                int64_t ln = stop - pos;
                if (clr) {
                    removed_any = 1;
                    if (open) {
                        open = 0;
                        frags->n++;
                    }
                    pending.n = 0;
                } else {
                    Frag *fr;
                    if (!open) {
                        fr = frag_open(frags);
                        fr->q0 = q;
                        fr->t0 = pos;
                        open = 1;
                    } else {
                        fr = &frags->a[frags->n];
                        for (int64_t g = 0; g < pending.n; g++) { // appended as they are, not merged
                            VEC_PUSH(fr->ops, pending.a[g]);
                        }
                    }
                    pending.n = 0;
                    push_op(&fr->ops, ln, op); // a matched op merges into the fragment's last op of its letter
                    add_query_iv(kept, r, q, ln);
                }
                q += ln;
                pos = stop;
            }
            tt += n;
            // (zero-length matched ops add nothing, as in the reference)
        } else if (op == 'I') {
            if (open) {
                Op o = { n, op };
                VEC_PUSH(pending, o);
            }
            q += n;
        } else {
            if (open) {
                Op o = { n, op };
                VEC_PUSH(pending, o);
            }
            tt += n;
        }
    }
    if (open) {
        frags->n++;
    }
    VEC_FREE(pending);
    if (!removed_any) {
        kept->n = kept0;
        return 0;
    }
    // the record's matched query intervals (all of them) are candidates for the query BED
    int64_t qq = 0;
    for (int64_t k = 0; k < ops->n; k++) {
        char op = ops->a[k].op;
        if (op == 'M' || op == '=' || op == 'X') {
            add_query_iv(removed, r, qq, ops->a[k].len);
            qq += ops->a[k].len;
        } else if (op == 'I') {
            qq += ops->a[k].len;
        }
    }
    return 1;
}

static void write_fragment(FILE *out, Fields *fl, Rec *r, Frag *fr, double m_ratio) {
    int64_t ql = 0, tl = 0, blk = 0, eq = 0, mm = 0;
    for (int64_t k = 0; k < fr->ops.n; k++) {
        Op *o = &fr->ops.a[k];
        if (o->op != 'D') ql += o->len;
        if (o->op != 'I') tl += o->len;
        blk += o->len;
        if (o->op == '=') eq += o->len;
        if (o->op == 'M') mm += o->len;
    }
    int64_t fqs, fqe;
    if (!r->rev) { fqs = r->qs + fr->q0; fqe = fqs + ql; } else { fqs = r->qe - fr->q0 - ql; fqe = r->qe - fr->q0; }
    int64_t matches = eq + (mm ? (int64_t)nearbyint((double)mm * m_ratio) : 0); // round half to even
    char **t = fl->f;
    fprintf(out, "%s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%s\t%s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%s",
            t[0], t[1], fqs, fqe, t[4], t[5], t[6], fr->t0, fr->t0 + tl, matches, blk, t[11]);
    for (int64_t k = 12; k < fl->n; k++) {
        if (k == r->cg_idx) {
            fputs("\tcg:Z:", out);
            for (int64_t g = 0; g < fr->ops.n; g++) {
                fprintf(out, "%" PRIi64 "%c", fr->ops.a[g].len, fr->ops.a[g].op);
            }
            continue;
        }
        const char *x = t[k];
        if (starts_with(x, "NM:") || starts_with(x, "AS:") || starts_with(x, "de:") || starts_with(x, "dv:") ||
            starts_with(x, "cs:") || starts_with(x, "ds:")) {
            continue;
        }
        fputc('\t', out);
        fputs(x, out);
    }
    fputc('\n', out);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Command line
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void usage(void) {
    fprintf(stderr, "paffy unanchor -i IN.paf -n NODES.fa -g GRAPH.gfa -r REF [options], version 0.1\n");
    fprintf(stderr, "Remove the hub (minigraph) anchors of every genome at tandem-repeat loci where the genomes' anchors disagree,\n"
                    "so that cactus realigns those loci. Input is one chromosome's PAF after filter_paf, every target a hub node\n");
    fprintf(stderr, " inputs\n");
    fprintf(stderr, "-i --inputFile FILE : Chromosome PAF (a file: it is read several times)\n");
    fprintf(stderr, "-n --nodes FILE : Hub node FASTA (headers sN or id=HUB|sN), plain or gzip\n");
    fprintf(stderr, "-g --gfa FILE : Minigraph rGFA, plain or gzip (S lines: LN, SN, SO, SR tags)\n");
    fprintf(stderr, "-r --ref NAME : Reference genome (PAF query prefix id=NAME|) [required]\n");
    fprintf(stderr, "-H --hub NAME : Hub genome (PAF target prefix id=NAME|) [_MINIGRAPH_]\n");
    fprintf(stderr, " outputs (each written only if given)\n");
    fprintf(stderr, "-o --outputFile FILE : Cut PAF [stdout]\n");
    fprintf(stderr, "-L --lociBed FILE : Cleared loci on the reference\n");
    fprintf(stderr, "-b --hubBed FILE : Cleared hub intervals\n");
    fprintf(stderr, "-v --vetoBed FILE : Candidate loci not cleared, with the reasons\n");
    fprintf(stderr, "-q --queryBed FILE : Query bases that lost every hub anchor (query name, start, end)\n");
    fprintf(stderr, "-c --callsBed FILE : Tandem calls on the reference path\n");
    fprintf(stderr, "-s --svReport FILE : SV features of the alt nodes of every bubble candidate locus\n");
    fprintf(stderr, "-e --logFile FILE : Summary (also on stderr at log level INFO)\n");
    fprintf(stderr, "-d --dryRun : Decide and write the BEDs, but copy the PAF unchanged (queryBed empty)\n");
    fprintf(stderr, " parameters\n");
    fprintf(stderr, "--gate comp|variable|none [comp]  --minHaps N [2]\n");
    fprintf(stderr, "--maxLen N : cactus-align's bandingLimit, sets the next two defaults [10000]\n");
    fprintf(stderr, "--maxSpan N [maxLen/2]  --maxString N [0.8*maxLen]  --cafTrim N [3]\n");
    fprintf(stderr, "--minPeriod N [2]  --maxPeriod N [100, at most 100]  --window N [20]  --minIdent F [0.8]\n");
    fprintf(stderr, "--pad N [10]  --padPeriods N [1]  --merge N [50]  --arrayMerge N [200]  --arrayMin N [10000]\n");
    fprintf(stderr, "--noSvVeto : Turn the SV veto off\n");
    fprintf(stderr, "--svMinBp N [50]  --svMaxLocal F [0.2]  --svInvMin F [0.5]  --svInvMargin F [0.3]  --svK N [13, at most 16]\n");
    fprintf(stderr, "--svTight N [100]  --svWide N [2000]\n");
    fprintf(stderr, "--keepCigarNoBoundary : Also clear cigar-class loci with no rank-0 node boundary (the v2 pilot)\n");
    fprintf(stderr, "--scanFasta FILE : Scan only: tandem-scan every sequence of FILE and write --callsBed\n");
    fprintf(stderr, "-t --threads N : Scan threads (output identical for any N) [1]\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

enum {
    O_GATE = 1000, O_MINHAPS, O_MAXLEN, O_MAXSPAN, O_MAXSTRING, O_CAFTRIM, O_MINPERIOD, O_MAXPERIOD, O_WINDOW,
    O_MINIDENT, O_PAD, O_PADPERIODS, O_MERGE, O_ARRAYMERGE, O_ARRAYMIN, O_NOSVVETO, O_SVMINBP, O_SVMAXLOCAL,
    O_SVINVMIN, O_SVINVMARGIN, O_SVK, O_SVTIGHT, O_SVWIDE, O_KEEPCNB, O_SCANFASTA
};

static double parse_double(const char *s, const char *what) {
    char *end;
    double v = strtod(s, &end);
    if (end == s || *end != '\0') {
        fprintf(stderr, "paffy unanchor: could not parse %s: %s\n", what, s);
        exit(1);
    }
    return v;
}

static int64_t parse_opt_int(const char *s, const char *what) {
    char *end;
    long long v = strtoll(s, &end, 10);
    if (end == s || *end != '\0') {
        fprintf(stderr, "paffy unanchor: could not parse %s: %s\n", what, s);
        exit(1);
    }
    return (int64_t)v;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Writing helpers
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static void write_log_file(Params *a) {
    if (a->log_file == NULL) return;
    FILE *f = open_out(a->log_file);
    for (int64_t i = 0; log_lines != NULL && i < stList_length(log_lines); i++) {
        fprintf(f, "%s\n", (char *)stList_get(log_lines, i));
    }
    close_out(f, a->log_file);
}

static void copy_file(const char *in_path, const char *out_path, int64_t *nrec) {
    FILE *in = fopen(in_path, "r");
    if (in == NULL) st_errnoAbort("Could not open %s", in_path);
    FILE *out = out_path == NULL ? stdout : open_out(out_path);
    char *buf = NULL;
    size_t cap = 0;
    ssize_t n;
    *nrec = 0;
    while ((n = getline(&buf, &cap, in)) >= 0) {
        if (is_blank_line(buf, n)) continue;
        (*nrec)++;
        fwrite(buf, 1, (size_t)n, out);
        if (n == 0 || buf[n - 1] != '\n') fputc('\n', out);
    }
    free(buf);
    fclose(in);
    if (out_path != NULL) close_out(out, out_path);
}

static const char *loci_header = "#contig\tstart\tend\tlocus\tperiod\tclass\tboundary\tnHaps\tcompHaps\tinsideDeps\t"
                                 "compAltnode\tcompReplace\tcompCigar\tminString\tmaxString\trefIntervals\taltNodes\t"
                                 "alleles\tcontigEnds\n";
static const char *hub_header = "#target\tstart\tend\tlocus\tkind\tclass\tboundary\tnHaps\tcompHaps\tmaxString\n";
static const char *sv_header = "#contig\tstart\tend\tlocus\tperiod\tstatus\tcompHaps\tnode\tlen\tfwdPasses\trevPasses\t"
                               "contigs\tgenomes\tkF\tkFw\tkRw\tnonTR\tsvReasons\n";

static void write_header_only(const char *path, const char *header) {
    if (path == NULL) return;
    FILE *f = open_out(path);
    fputs(header, f);
    close_out(f, path);
}

static const char *veto_header(Params *a) {
    return a->sv_veto ? "#contig\tstart\tend\tperiod\treason\tallReasons\tcompHaps\tsvDetail\n"
                      : "#contig\tstart\tend\tperiod\treason\tallReasons\tcompHaps\n";
}

static void to_local(State *st, int64_t x, const char **contig, int64_t *local) {
    int64_t i = upper_bound(st->coff, st->ncontigs, x) - 1;
    *contig = st->contigs[i];
    *local = x - st->coff[i];
}

static void reasons_string(uint32_t mask, CharVec *out) {
    out->n = 0;
    cv_str(out);
    bool first = 1;
    for (int r = 0; r < R_NUM; r++) {
        if (mask & (1u << r)) {
            cv_printf(out, "%s%s", first ? "" : ",", reason_names[r]);
            first = 0;
        }
    }
    if (first) cv_printf(out, "-");
}

static int first_reason(uint32_t mask) {
    for (int r = 0; r < R_NUM; r++) {
        if (mask & (1u << r)) return r;
    }
    return -1;
}

typedef struct { int64_t base, rank, node; } GNode;

static int cmp_gnode(const void *x, const void *y) {
    const GNode *a = x, *b = y;
    if (a->base != b->base) return a->base < b->base ? -1 : 1;
    return a->rank < b->rank ? -1 : (a->rank > b->rank ? 1 : 0);
}

typedef struct { int64_t base, order, node; } SeqNode;

static int cmp_seqnode(const void *x, const void *y) {
    const SeqNode *a = x, *b = y;
    if (a->base != b->base) return a->base < b->base ? -1 : 1;
    return a->order < b->order ? -1 : (a->order > b->order ? 1 : 0);
}

static int cmp_pass(const void *x, const void *y) {
    const Pass *a = x, *b = y;
    if (a->li != b->li) return a->li < b->li ? -1 : 1;
    return a->order < b->order ? -1 : (a->order > b->order ? 1 : 0);
}

static NodeTable *sort_nt; // for the node-rank comparators

static int cmp_altin(const void *x, const void *y) {
    const AltIn *a = x, *b = y;
    int64_t ra = sort_nt->v.a[a->node].rank, rb = sort_nt->v.a[b->node].rank;
    if (ra != rb) return ra < rb ? -1 : 1;
    return a->li < b->li ? -1 : (a->li > b->li ? 1 : 0);
}

static int cmp_hubiv_locus(const void *x, const void *y) {
    const HubIv *a = x, *b = y;
    if (a->li != b->li) return a->li < b->li ? -1 : 1;
    int64_t ra = sort_nt->v.a[a->node].rank, rb = sort_nt->v.a[b->node].rank;
    if (ra != rb) return ra < rb ? -1 : 1;
    if (a->ts != b->ts) return a->ts < b->ts ? -1 : 1;
    if (a->te != b->te) return a->te < b->te ? -1 : 1;
    return (int)a->alt - (int)b->alt; // 'alt' < 'ref'
}

static int cmp_hubiv_target(const void *x, const void *y) {
    const HubIv *a = x, *b = y;
    int64_t ra = sort_nt->v.a[a->node].rank, rb = sort_nt->v.a[b->node].rank;
    if (ra != rb) return ra < rb ? -1 : 1;
    if (a->ts != b->ts) return a->ts < b->ts ? -1 : 1;
    if (a->te != b->te) return a->te < b->te ? -1 : 1;
    if (a->li != b->li) return a->li < b->li ? -1 : 1;
    return (int)a->alt - (int)b->alt;
}

static int cmp_name_index(const void *x, const void *y) {
    int64_t a = *(const int64_t *)x, b = *(const int64_t *)y;
    return strcmp(sort_nt->v.a[a].name, sort_nt->v.a[b].name);
}

typedef struct { int64_t k, contig; } KC;

static int cmp_kc(const void *x, const void *y) {
    const KC *a = x, *b = y;
    if (a->k != b->k) return a->k < b->k ? -1 : 1;
    return a->contig < b->contig ? -1 : (a->contig > b->contig ? 1 : 0);
}

static char **qsort_names; // contig names for cmp_qc_index

static int cmp_qc_index(const void *x, const void *y) {
    return strcmp(qsort_names[*(const int64_t *)x], qsort_names[*(const int64_t *)y]);
}

typedef VEC(char *) StrVec;

// free what the inputs were read into (shared by the normal end and the nothing-to-clear exit)
static void free_inputs(State *st, Scratch *sc, StrVec *ref_names, I64Vec *ref_lens, stHash *ref_len_index) {
    for (int64_t i = 0; i < ref_names->n; i++) free(ref_names->a[i]);
    VEC_FREE(*ref_names);
    VEC_FREE(*ref_lens);
    stHash_destruct(ref_len_index);
    for (int64_t i = 0; i < st->qc.n; i++) {
        free(st->qc.a[i].name);
        VEC_FREE(st->qc.a[i].recs);
    }
    VEC_FREE(st->qc);
    stHash_destruct(st->qindex);
    stHash_destruct(st->gindex);
    for (int64_t i = 0; i < st->nt.v.n; i++) {
        free(st->nt.v.a[i].name);
        free(st->nt.v.a[i].rank0_contig);
        free(st->nt.v.a[i].seq);
    }
    VEC_FREE(st->nt.v);
    stHash_destruct(st->nt.index);
    VEC_FREE(sc->line);
    free(sc->fl.f);
    VEC_FREE(sc->ops);
    VEC_FREE(sc->runs);
    free(st->ref_q);
    free(st->hub_t);
    if (log_lines != NULL) {
        stList_destruct(log_lines);
        log_lines = NULL;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Main
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int paffy_unanchor_main(int argc, char *argv[]) {
    double T0 = now_seconds();
    char *logLevelString = NULL;
    Params A;
    memset(&A, 0, sizeof(A));
    Params *a = &A;
    a->hub = "_MINIGRAPH_";
    a->sv_veto = 1;
    a->gate = GATE_COMP;
    a->min_haps = 2;
    a->max_len = 10000;
    a->max_span = -1;
    a->max_string = -1;
    a->caf_trim = 3;
    a->min_period = 2;
    a->max_period = 100;
    a->window = 20;
    a->min_ident = 0.8;
    a->pad = 10;
    a->pad_periods = 1;
    a->merge = 50;
    a->array_merge = 200;
    a->array_min = 10000;
    a->sv_min_bp = 50;
    a->sv_max_local = 0.2;
    a->sv_inv_min = 0.5;
    a->sv_inv_margin = 0.3;
    a->sv_k = 13;
    a->sv_tight = 100;
    a->sv_wide = 2000;
    a->threads = 1;

    while (1) {
        static struct option long_options[] = {
            { "logLevel", required_argument, 0, 'l' }, { "inputFile", required_argument, 0, 'i' },
            { "nodes", required_argument, 0, 'n' }, { "gfa", required_argument, 0, 'g' },
            { "ref", required_argument, 0, 'r' }, { "hub", required_argument, 0, 'H' },
            { "outputFile", required_argument, 0, 'o' }, { "lociBed", required_argument, 0, 'L' },
            { "hubBed", required_argument, 0, 'b' }, { "vetoBed", required_argument, 0, 'v' },
            { "queryBed", required_argument, 0, 'q' }, { "callsBed", required_argument, 0, 'c' },
            { "svReport", required_argument, 0, 's' }, { "logFile", required_argument, 0, 'e' },
            { "dryRun", no_argument, 0, 'd' }, { "threads", required_argument, 0, 't' },
            { "help", no_argument, 0, 'h' },
            { "gate", required_argument, 0, O_GATE }, { "minHaps", required_argument, 0, O_MINHAPS },
            { "maxLen", required_argument, 0, O_MAXLEN }, { "maxSpan", required_argument, 0, O_MAXSPAN },
            { "maxString", required_argument, 0, O_MAXSTRING }, { "cafTrim", required_argument, 0, O_CAFTRIM },
            { "minPeriod", required_argument, 0, O_MINPERIOD }, { "maxPeriod", required_argument, 0, O_MAXPERIOD },
            { "window", required_argument, 0, O_WINDOW }, { "minIdent", required_argument, 0, O_MINIDENT },
            { "pad", required_argument, 0, O_PAD }, { "padPeriods", required_argument, 0, O_PADPERIODS },
            { "merge", required_argument, 0, O_MERGE }, { "arrayMerge", required_argument, 0, O_ARRAYMERGE },
            { "arrayMin", required_argument, 0, O_ARRAYMIN }, { "noSvVeto", no_argument, 0, O_NOSVVETO },
            { "svMinBp", required_argument, 0, O_SVMINBP }, { "svMaxLocal", required_argument, 0, O_SVMAXLOCAL },
            { "svInvMin", required_argument, 0, O_SVINVMIN }, { "svInvMargin", required_argument, 0, O_SVINVMARGIN },
            { "svK", required_argument, 0, O_SVK }, { "svTight", required_argument, 0, O_SVTIGHT },
            { "svWide", required_argument, 0, O_SVWIDE }, { "keepCigarNoBoundary", no_argument, 0, O_KEEPCNB },
            { "scanFasta", required_argument, 0, O_SCANFASTA },
            { 0, 0, 0, 0 } };
        int option_index = 0;
        int key = getopt_long(argc, argv, "l:i:n:g:r:H:o:L:b:v:q:c:s:e:dt:h", long_options, &option_index);
        if (key == -1) break;
        switch (key) {
            case 'l': logLevelString = optarg; break;
            case 'i': a->paf = optarg; break;
            case 'n': a->nodes = optarg; break;
            case 'g': a->gfa = optarg; break;
            case 'r': a->ref = optarg; break;
            case 'H': a->hub = optarg; break;
            case 'o': a->out_paf = optarg; break;
            case 'L': a->loci_bed = optarg; break;
            case 'b': a->hub_bed = optarg; break;
            case 'v': a->veto_bed = optarg; break;
            case 'q': a->query_bed = optarg; break;
            case 'c': a->calls_bed = optarg; break;
            case 's': a->sv_report = optarg; break;
            case 'e': a->log_file = optarg; break;
            case 'd': a->dry_run = 1; break;
            case 't': a->threads = parse_opt_int(optarg, "--threads"); break;
            case 'h': usage(); return 0;
            case O_GATE:
                if (strcmp(optarg, "comp") == 0) a->gate = GATE_COMP;
                else if (strcmp(optarg, "variable") == 0) a->gate = GATE_VARIABLE;
                else if (strcmp(optarg, "none") == 0) a->gate = GATE_NONE;
                else { fprintf(stderr, "paffy unanchor: --gate must be comp, variable or none\n"); return 1; }
                break;
            case O_MINHAPS: a->min_haps = parse_opt_int(optarg, "--minHaps"); break;
            case O_MAXLEN: a->max_len = parse_opt_int(optarg, "--maxLen"); break;
            case O_MAXSPAN: a->max_span = parse_opt_int(optarg, "--maxSpan"); break;
            case O_MAXSTRING: a->max_string = parse_opt_int(optarg, "--maxString"); break;
            case O_CAFTRIM: a->caf_trim = parse_opt_int(optarg, "--cafTrim"); break;
            case O_MINPERIOD: a->min_period = parse_opt_int(optarg, "--minPeriod"); break;
            case O_MAXPERIOD: a->max_period = parse_opt_int(optarg, "--maxPeriod"); break;
            case O_WINDOW: a->window = parse_opt_int(optarg, "--window"); break;
            case O_MINIDENT: a->min_ident = parse_double(optarg, "--minIdent"); break;
            case O_PAD: a->pad = parse_opt_int(optarg, "--pad"); break;
            case O_PADPERIODS: a->pad_periods = parse_opt_int(optarg, "--padPeriods"); break;
            case O_MERGE: a->merge = parse_opt_int(optarg, "--merge"); break;
            case O_ARRAYMERGE: a->array_merge = parse_opt_int(optarg, "--arrayMerge"); break;
            case O_ARRAYMIN: a->array_min = parse_opt_int(optarg, "--arrayMin"); break;
            case O_NOSVVETO: a->sv_veto = 0; break;
            case O_SVMINBP: a->sv_min_bp = parse_opt_int(optarg, "--svMinBp"); break;
            case O_SVMAXLOCAL: a->sv_max_local = parse_double(optarg, "--svMaxLocal"); break;
            case O_SVINVMIN: a->sv_inv_min = parse_double(optarg, "--svInvMin"); break;
            case O_SVINVMARGIN: a->sv_inv_margin = parse_double(optarg, "--svInvMargin"); break;
            case O_SVK: a->sv_k = parse_opt_int(optarg, "--svK"); break;
            case O_SVTIGHT: a->sv_tight = parse_opt_int(optarg, "--svTight"); break;
            case O_SVWIDE: a->sv_wide = parse_opt_int(optarg, "--svWide"); break;
            case O_KEEPCNB: a->keep_cigar_no_boundary = 1; break;
            case O_SCANFASTA: a->scan_fasta = optarg; break;
            default: usage(); return 1;
        }
    }
    if (optind < argc) {
        fprintf(stderr, "paffy unanchor: unexpected argument %s\n", argv[optind]);
        usage();
        return 1;
    }
    st_setLogLevelFromString(logLevelString);
    if (a->max_span < 0) a->max_span = a->max_len / 2;
    if (a->max_string < 0) a->max_string = (int64_t)(0.8 * (double)a->max_len);
    if (a->max_period > 100 || a->min_period < 1 || a->min_period > a->max_period) {
        fprintf(stderr, "paffy unanchor: periods must satisfy 1 <= minPeriod <= maxPeriod <= 100\n");
        return 1;
    }
    if (a->sv_k < 1 || a->sv_k > 16) {
        fprintf(stderr, "paffy unanchor: --svK must be between 1 and 16\n");
        return 1;
    }
    if (a->threads < 1) a->threads = 1;
    if (a->scan_fasta != NULL) {
        if (a->calls_bed == NULL) {
            fprintf(stderr, "paffy unanchor: --scanFasta needs --callsBed\n");
            return 1;
        }
        int r = scan_only(a);
        write_log_file(a);
        if (log_lines != NULL) {
            stList_destruct(log_lines);
            log_lines = NULL;
        }
        return r;
    }
    if (a->paf == NULL || a->nodes == NULL || a->gfa == NULL || a->ref == NULL) {
        fprintf(stderr, "paffy unanchor: --inputFile, --nodes, --gfa and --ref are required\n");
        usage();
        return 1;
    }

    {
        CharVec cl = { 0 };
        cv_printf(&cl, "paffy unanchor");
        for (int i = 1; i < argc; i++) cv_printf(&cl, " %s", argv[i]);
        ulog("%s", cv_str(&cl));
        VEC_FREE(cl);
    }
    ulog("params arrayMerge=%" PRIi64 " arrayMin=%" PRIi64 " cafTrim=%" PRIi64 " dryRun=%d gate=%s hub=%s "
         "keepCigarNoBoundary=%d maxLen=%" PRIi64 " maxPeriod=%" PRIi64 " maxSpan=%" PRIi64 " maxString=%" PRIi64
         " merge=%" PRIi64 " minHaps=%" PRIi64 " minIdent=%g minPeriod=%" PRIi64 " pad=%" PRIi64 " padPeriods=%"
         PRIi64 " ref=%s svInvMargin=%g svInvMin=%g svK=%" PRIi64 " svMaxLocal=%g svMinBp=%" PRIi64 " svReport=%d "
         "svTight=%" PRIi64 " svVeto=%d svWide=%" PRIi64 " threads=%" PRIi64 " window=%" PRIi64,
         a->array_merge, a->array_min, a->caf_trim, (int)a->dry_run,
         a->gate == GATE_COMP ? "comp" : a->gate == GATE_VARIABLE ? "variable" : "none", a->hub,
         (int)a->keep_cigar_no_boundary, a->max_len, a->max_period, a->max_span, a->max_string, a->merge,
         a->min_haps, a->min_ident, a->min_period, a->pad, a->pad_periods, a->ref, a->sv_inv_margin, a->sv_inv_min,
         a->sv_k, a->sv_max_local, a->sv_min_bp, a->sv_report != NULL, a->sv_tight, (int)a->sv_veto, a->sv_wide,
         a->threads, a->window);

    State S;
    memset(&S, 0, sizeof(S));
    State *st = &S;
    st->a = a;
    st->ref_q = strprintf("id=%s|", a->ref);
    st->hub_t = strprintf("id=%s|", a->hub);
    st->nt.index = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    st->qindex = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    st->gindex = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    size_t ref_q_len = strlen(st->ref_q), hub_t_len = strlen(st->hub_t);

    //////////////////////////////////////////////
    // (0) GFA S lines and the node FASTA
    //////////////////////////////////////////////

    int64_t n_s, n_rank0, n_alt;
    read_gfa(st, &n_s, &n_rank0, &n_alt);
    int64_t n_fasta = 0;
    {
        FaVec fa = { 0 };
        read_fasta(a->nodes, &fa);
        for (int64_t i = 0; i < fa.n; i++) {
            char *name = starts_with(fa.a[i].name, st->hub_t) ? stString_copy(fa.a[i].name)
                                                               : strprintf("%s%s", st->hub_t, fa.a[i].name);
            int64_t j = node_get(&st->nt, name, 1);
            free(name);
            Node *nd = &st->nt.v.a[j];
            if (!(nd->flags & NODE_FASTA)) n_fasta++;
            free(nd->seq);
            nd->seq = fa.a[i].seq;
            nd->seq_len = fa.a[i].len;
            nd->flags |= NODE_FASTA;
            free(fa.a[i].name);
        }
        VEC_FREE(fa);
    }

    //////////////////////////////////////////////
    // (1) PAF pass 1: check, target lengths, the reference's runs, the query contigs' records
    //////////////////////////////////////////////

    FILE *in = fopen(a->paf, "r");
    if (in == NULL) {
        st_errnoAbort("Could not open the input PAF %s", a->paf);
    }
    Scratch sc;
    memset(&sc, 0, sizeof(sc));
    stHash *ref_len_index = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL, NULL);
    StrVec ref_names = { 0 };
    I64Vec ref_lens = { 0 };
    typedef struct { int64_t contig, g, n, node, t0; bool rev; } RefRun;
    VEC(RefRun) ref_runs = { 0 };
    int64_t nrec = 0, bad_target = 0;
    {
        char *buf = NULL;
        size_t cap = 0;
        ssize_t n;
        off_t off = 0;
        while ((n = getline(&buf, &cap, in)) >= 0) {
            off_t line_off = off;
            off += n;
            if (is_blank_line(buf, n)) continue;
            if (buf[n - 1] == '\n') buf[--n] = '\0';
            int64_t ri = nrec++;
            VEC_GROW(sc.line, n + 1);
            memcpy(sc.line.a, buf, (size_t)n + 1);
            split_fields(sc.line.a, &sc.fl);
            Rec r;
            if (sc.fl.n < 12) {
                die2("PAF record with %" PRIi64 " fields (at least 12 needed): %s", sc.fl.n, sc.fl.f[0]);
            }
            if (!starts_with(sc.fl.f[5], st->hub_t) || strcmp(sc.fl.f[0], sc.fl.f[5]) == 0) {
                bad_target++;
                if (bad_target <= 5) {
                    ulog("ERROR: record target is not the hub (%s): %s %s", st->hub_t, sc.fl.f[0], sc.fl.f[5]);
                }
                continue;
            }
            parse_record(&sc.fl, &r, &sc.ops);
            int64_t node = node_get(&st->nt, r.t, 1);
            st->nt.v.a[node].tlen = r.tlen;
            if (starts_with(r.q, st->ref_q)) {
                const char *c = r.q + ref_q_len;
                void *x = stHash_search(ref_len_index, (void *)c);
                int64_t ci;
                if (x == NULL) {
                    char *cc = stString_copy(c);
                    VEC_PUSH(ref_names, cc);
                    VEC_PUSH(ref_lens, 0);
                    ci = ref_names.n - 1;
                    stHash_insert(ref_len_index, cc, (void *)(intptr_t)(ci + 1));
                } else {
                    ci = (int64_t)(intptr_t)x - 1;
                }
                ref_lens.a[ci] = r.qlen;
                if (!r.secondary) {
                    record_runs(&r, &sc.ops, &sc.runs);
                    for (int64_t k = 0; k < sc.runs.n; k++) {
                        RefRun rr = { ci, sc.runs.a[k].qf, sc.runs.a[k].n, node, sc.runs.a[k].t0, r.rev };
                        VEC_PUSH(ref_runs, rr);
                    }
                }
            } else if (r.secondary) {
                // secondaries do not walk a haplotype; an alt node one lands on counts as used outside any locus
                if (!(st->nt.v.a[node].flags & NODE_RANK0)) st->nt.v.a[node].flags |= NODE_ALT_OUT;
            } else {
                void *x = stHash_search(st->qindex, r.q);
                int64_t qi;
                if (x == NULL) {
                    QContig qc;
                    memset(&qc, 0, sizeof(qc));
                    qc.name = stString_copy(r.q);
                    VEC_PUSH(st->qc, qc);
                    qi = st->qc.n - 1;
                    stHash_insert(st->qindex, st->qc.a[qi].name, (void *)(intptr_t)(qi + 1));
                } else {
                    qi = (int64_t)(intptr_t)x - 1;
                }
                VEC_PUSH(st->qc.a[qi].recs, (int64_t)line_off);
                VEC_PUSH(st->qc.a[qi].recs, ri);
            }
        }
        free(buf);
        if (ferror(in)) st_errnoAbort("Error reading %s", a->paf);
    }
    if (bad_target) {
        ulog("ERROR: %" PRIi64 " records do not target the hub genome %s (self-alignments from <graphmap collapse>, or a "
             "PAF from removeMinigraphFromPAF): refusing to run, since clearing keyed on hub positions cannot remove "
             "their anchors", bad_target, a->hub);
        write_log_file(a);
        exit(2);
    }
    ulog("PAF records %" PRIi64, nrec);

    // node names in byte order
    {
        int64_t nn = st->nt.v.n;
        int64_t *idx = st_malloc(sizeof(int64_t) * (size_t)(nn > 0 ? nn : 1));
        for (int64_t i = 0; i < nn; i++) idx[i] = i;
        sort_nt = &st->nt;
        qsort(idx, (size_t)nn, sizeof(int64_t), cmp_name_index);
        for (int64_t i = 0; i < nn; i++) st->nt.v.a[idx[i]].rank = i;
        free(idx);
    }
    int64_t nn = st->nt.v.n;
    Node *N = st->nt.v.a;
    ulog("GFA segments %" PRIi64 ": rank-0 of %s %" PRIi64 ", other %" PRIi64, n_s, a->ref, n_rank0, n_alt);
    if (n_fasta > 0) {
        int64_t miss = 0, badlen = 0, extra = 0;
        char *bad_example = NULL;
        for (int64_t i = 0; i < nn; i++) {
            bool in_gfa = N[i].flags & (NODE_RANK0 | NODE_GFA_ALT);
            if (N[i].flags & NODE_RANK0) {
                if (!(N[i].flags & NODE_FASTA)) miss++;
            }
            if (N[i].flags & NODE_GFA_ALT) {
                if (!(N[i].flags & NODE_FASTA)) miss++;
            }
            if (in_gfa && (N[i].flags & NODE_FASTA) && N[i].seq_len != N[i].ln) {
                badlen += ((N[i].flags & NODE_RANK0) != 0) + ((N[i].flags & NODE_GFA_ALT) != 0);
                if (bad_example == NULL) bad_example = N[i].name;
            }
            if ((N[i].flags & NODE_FASTA) && !in_gfa) extra++;
        }
        ulog("node FASTA vs GFA: %" PRIi64 " sequences; GFA segments missing from FASTA %" PRIi64 "; length mismatches %"
             PRIi64 "; FASTA-only %" PRIi64, n_fasta, miss, badlen, extra);
        if (badlen) {
            ulog("ERROR: node FASTA and GFA disagree on node lengths, e.g. %s", bad_example);
            write_log_file(a);
            exit(2);
        }
    }

    //////////////////////////////////////////////
    // (2) The reference layout and frame
    //////////////////////////////////////////////

    for (int64_t i = 0; i < nn; i++) {
        if (!(N[i].flags & NODE_RANK0)) continue;
        const char *c = N[i].rank0_contig;
        void *x = stHash_search(ref_len_index, (void *)c);
        int64_t ci;
        if (x == NULL) {
            char *cc = stString_copy(c);
            VEC_PUSH(ref_names, cc);
            VEC_PUSH(ref_lens, 0);
            ci = ref_names.n - 1;
            stHash_insert(ref_len_index, cc, (void *)(intptr_t)(ci + 1));
        } else {
            ci = (int64_t)(intptr_t)x - 1;
        }
        // the contig is at least as long as its rank-0 nodes reach (max with the PAF query length)
        int64_t reach = N[i].so + N[i].ln;
        if (reach > ref_lens.a[ci]) ref_lens.a[ci] = reach;
    }
    st->ncontigs = ref_names.n;
    if (st->ncontigs == 0) {
        ulog("WARNING: no reference contig found (no records with query prefix %s and no rank-0 GFA nodes): nothing "
             "to clear, the PAF is copied unchanged", st->ref_q);
        int64_t nout;
        copy_file(a->paf, a->out_paf, &nout);
        write_header_only(a->loci_bed, loci_header);
        write_header_only(a->hub_bed, hub_header);
        write_header_only(a->veto_bed, veto_header(a));
        write_header_only(a->query_bed, "");
        write_header_only(a->calls_bed, "");
        if (a->sv_report != NULL) write_header_only(a->sv_report, sv_header);
        write_log_file(a);
        fclose(in);
        VEC_FREE(ref_runs);
        free_inputs(st, &sc, &ref_names, &ref_lens, ref_len_index);
        return 0;
    }
    // contigs in byte order
    int64_t *corder = st_malloc(sizeof(int64_t) * (size_t)st->ncontigs);
    for (int64_t i = 0; i < st->ncontigs; i++) corder[i] = i;
    qsort_names = ref_names.a;
    qsort(corder, (size_t)st->ncontigs, sizeof(int64_t), cmp_qc_index);
    int64_t *cpos = st_malloc(sizeof(int64_t) * (size_t)st->ncontigs); // original index -> sorted position
    st->contigs = st_malloc(sizeof(char *) * (size_t)st->ncontigs);
    st->clen = st_malloc(sizeof(int64_t) * (size_t)st->ncontigs);
    st->coff = st_malloc(sizeof(int64_t) * (size_t)st->ncontigs);
    int64_t SEP = 2 * (a->max_period + (a->window > a->max_period ? a->window : a->max_period)) + 1;
    {
        int64_t tot = 0;
        for (int64_t k = 0; k < st->ncontigs; k++) {
            int64_t i = corder[k];
            cpos[i] = k;
            st->contigs[k] = ref_names.a[i];
            st->clen[k] = ref_lens.a[i];
            st->coff[k] = tot;
            tot += ref_lens.a[i] + SEP;
        }
        st->G = tot;
    }
    {
        CharVec m = { 0 };
        cv_printf(&m, "reference contigs: ");
        for (int64_t k = 0; k < st->ncontigs; k++) {
            cv_printf(&m, "%s%s (%" PRIi64 " bp)", k ? ", " : "", st->contigs[k], st->clen[k]);
        }
        ulog("%s", cv_str(&m));
        VEC_FREE(m);
    }
    // node bases, the reference-path sequence and the node boundaries
    st->seq = st_malloc((size_t)st->G + 1);
    memset(st->seq, 'N', (size_t)st->G);
    st->seq[st->G] = '\0';
    int64_t missing_bases = 0, framed = 0;
    {
        VEC(SeqNode) sn = { 0 };
        for (int64_t i = 0; i < nn; i++) {
            if (!(N[i].flags & NODE_RANK0)) continue;
            void *x = stHash_search(ref_len_index, N[i].rank0_contig);
            int64_t k = cpos[(int64_t)(intptr_t)x - 1];
            N[i].base = st->coff[k] + N[i].so;
            framed++;
            SeqNode s = { N[i].base, N[i].gfa_order, i };
            VEC_PUSH(sn, s);
        }
        if (sn.n > 0) qsort(sn.a, (size_t)sn.n, sizeof(SeqNode), cmp_seqnode);
        for (int64_t k = 0; k < sn.n; k++) {
            Node *nd = &N[sn.a[k].node];
            if (!(nd->flags & NODE_FASTA)) {
                missing_bases += nd->ln;
            } else {
                if (nd->base < 0 || nd->base + nd->seq_len > st->G) {
                    die2("rank-0 node %s lies outside its reference contig", nd->name);
                }
                memcpy(st->seq + nd->base, nd->seq, (size_t)nd->seq_len);
            }
            if (nd->so > 0) VEC_PUSH(st->bounds, nd->base);
        }
        VEC_FREE(sn);
        if (st->bounds.n > 0) qsort(st->bounds.a, (size_t)st->bounds.n, sizeof(int64_t), cmp_i64);
    }
    // confirmed positions: the reference's own primary runs place node offset t0 at its GFA position
    int64_t n_ref_alt = 0;
    IvVec conf = { 0 };
    for (int64_t k = 0; k < ref_runs.n; k++) {
        RefRun *rr = &ref_runs.a[k];
        int64_t g = st->coff[cpos[rr->contig]] + rr->g;
        Node *nd = &N[rr->node];
        if (nd->base < 0) {
            if (!(nd->flags & NODE_REF_ALT)) n_ref_alt++;
            nd->flags |= NODE_REF_ALT | NODE_ALT_OUT;
            continue;
        }
        if (!rr->rev && nd->base + rr->t0 == g) {
            Iv x = { g, g + rr->n };
            VEC_PUSH(conf, x);
        }
    }
    VEC_FREE(ref_runs);
    if (conf.n > 0) qsort(conf.a, (size_t)conf.n, sizeof(Iv), cmp_iv);
    IvVec cov2 = { 0 }, confm = { 0 }, unframed = { 0 };
    for (int64_t k = 0; k < conf.n; k++) {
        int64_t s = conf.a[k].s, e = conf.a[k].e;
        if (confm.n > 0 && s < confm.a[confm.n - 1].e) {
            Iv x = { s, e < confm.a[confm.n - 1].e ? e : confm.a[confm.n - 1].e };
            VEC_PUSH(cov2, x);
            if (e > confm.a[confm.n - 1].e) confm.a[confm.n - 1].e = e;
        } else if (confm.n > 0 && s == confm.a[confm.n - 1].e) {
            confm.a[confm.n - 1].e = e;
        } else {
            Iv x = { s, e };
            VEC_PUSH(confm, x);
        }
    }
    VEC_FREE(conf);
    for (int64_t k = 0; k < st->ncontigs; k++) {
        int64_t lo = st->coff[k], hi = st->coff[k] + st->clen[k], cur = lo;
        for (int64_t j = 0; j < confm.n; j++) {
            if (confm.a[j].e <= lo || confm.a[j].s >= hi) continue;
            if (confm.a[j].s > cur) {
                Iv x = { cur, confm.a[j].s };
                VEC_PUSH(unframed, x);
            }
            if (confm.a[j].e > cur) cur = confm.a[j].e;
        }
        if (cur < hi) {
            Iv x = { cur, hi };
            VEC_PUSH(unframed, x);
        }
    }
    int64_t cov2_bp = 0;
    for (int64_t k = 0; k < cov2.n; k++) {
        VEC_PUSH(unframed, cov2.a[k]);
        cov2_bp += cov2.a[k].e - cov2.a[k].s;
    }
    VEC_FREE(cov2);
    VEC_FREE(confm);
    if (unframed.n > 0) qsort(unframed.a, (size_t)unframed.n, sizeof(Iv), cmp_iv);
    int64_t unframed_bp = 0;
    {
        int64_t m = -1;
        for (int64_t k = 0; k < unframed.n; k++) {
            unframed_bp += unframed.a[k].e - unframed.a[k].s;
            if (unframed.a[k].e > m) m = unframed.a[k].e;
            VEC_PUSH(st->uf_s, unframed.a[k].s);
            VEC_PUSH(st->uf_emax, m);
        }
    }
    ulog("frame (gfa): framed nodes %" PRIi64 "; reference-path bases missing from the node FASTA %" PRIi64
         "; unframed reference bp %" PRIi64 " in %" PRIi64 " intervals (doubly placed %" PRIi64 " bp); reference runs "
         "on non-rank-0 nodes %" PRIi64 " nodes", framed, missing_bases, unframed_bp, unframed.n, cov2_bp, n_ref_alt);
    VEC_FREE(unframed);

    //////////////////////////////////////////////
    // Tandem scan, one contig at a time
    //////////////////////////////////////////////

    double t1 = now_seconds();
    I64Vec calls = { 0 }; // (start, end, period) triples, global
    for (int64_t k = 0; k < st->ncontigs; k++) {
        int64_t nc;
        int64_t *c = tandem_scan(st->seq + st->coff[k], st->clen[k], a->min_period, a->max_period, a->window,
                                 a->min_ident, a->threads, &nc);
        for (int64_t j = 0; j < nc; j++) {
            VEC_PUSH(calls, c[3 * j] + st->coff[k]);
            VEC_PUSH(calls, c[3 * j + 1] + st->coff[k]);
            VEC_PUSH(calls, c[3 * j + 2]);
        }
        free(c);
    }
    int64_t ncalls = calls.n / 3, calls_bp = 0;
    for (int64_t j = 0; j < ncalls; j++) calls_bp += calls.a[3 * j + 1] - calls.a[3 * j];
    // contigs are laid out in order and each contig's calls are sorted, so the global list is sorted
    ulog("tandem scan periods %" PRIi64 "-%" PRIi64 ": %" PRIi64 " calls, %" PRIi64 " bp, %.1fs", a->min_period,
         a->max_period, ncalls, calls_bp, now_seconds() - t1);
    if (a->calls_bed != NULL) {
        FILE *f = open_out(a->calls_bed);
        for (int64_t j = 0; j < ncalls; j++) {
            const char *c;
            int64_t ls;
            to_local(st, calls.a[3 * j], &c, &ls);
            fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\n", c, ls, ls + calls.a[3 * j + 1] - calls.a[3 * j],
                    calls.a[3 * j + 2]);
        }
        close_out(f, a->calls_bed);
    }

    //////////////////////////////////////////////
    // Arrays and candidate loci
    //////////////////////////////////////////////

    IvVec arrays = { 0 };
    for (int64_t j = 0; j < ncalls; j++) {
        int64_t s = calls.a[3 * j], e = calls.a[3 * j + 1];
        if (arrays.n > 0 && s - arrays.a[arrays.n - 1].e <= a->array_merge) {
            if (e > arrays.a[arrays.n - 1].e) arrays.a[arrays.n - 1].e = e;
        } else {
            Iv x = { s, e };
            VEC_PUSH(arrays, x);
        }
    }
    {
        int64_t j = 0;
        for (int64_t k = 0; k < arrays.n; k++) {
            if (arrays.a[k].e - arrays.a[k].s >= a->array_min) arrays.a[j++] = arrays.a[k];
        }
        arrays.n = j;
    }
    I64Vec astarts = { 0 };
    int64_t arrays_bp = 0;
    for (int64_t k = 0; k < arrays.n; k++) {
        VEC_PUSH(astarts, arrays.a[k].s);
        arrays_bp += arrays.a[k].e - arrays.a[k].s;
    }
    for (int64_t j = 0; j < ncalls; j++) {
        int64_t s = calls.a[3 * j], e = calls.a[3 * j + 1], p = calls.a[3 * j + 2];
        if (in_array(&astarts, &arrays, s, e)) continue;
        int64_t pad = a->pad + a->pad_periods * p;
        int64_t s2 = s - pad > 0 ? s - pad : 0, e2 = e + pad;
        int64_t ci = upper_bound(st->coff, st->ncontigs, s) - 1;
        if (s2 < st->coff[ci]) s2 = st->coff[ci];
        if (e2 > st->coff[ci] + st->clen[ci]) e2 = st->coff[ci] + st->clen[ci];
        if (st->loci.n > 0 && s2 - st->loci.a[st->loci.n - 1].e <= a->merge) {
            Locus *l = &st->loci.a[st->loci.n - 1];
            if (e2 > l->e) l->e = e2;
            if (p > l->p) l->p = p;
        } else {
            Locus l = { s2, e2, p };
            VEC_PUSH(st->loci, l);
        }
    }
    {
        int64_t j = 0;
        for (int64_t k = 0; k < st->loci.n; k++) {
            if (!in_array(&astarts, &arrays, st->loci.a[k].s, st->loci.a[k].e)) st->loci.a[j++] = st->loci.a[k];
        }
        st->loci.n = j;
    }
    int64_t narrays = arrays.n;
    VEC_FREE(calls);
    VEC_FREE(arrays);
    VEC_FREE(astarts);
    int64_t nl = st->loci.n, loci_bp = 0;
    Locus *L = st->loci.a;
    for (int64_t k = 0; k < nl; k++) {
        VEC_PUSH(st->lstarts, L[k].s);
        loci_bp += L[k].e - L[k].s;
    }
    ulog("arrays (calls merged within %" PRIi64 " bp, >= %" PRIi64 " bp): %" PRIi64 ", %" PRIi64 " bp", a->array_merge,
         a->array_min, narrays, arrays_bp);
    ulog("candidate loci %" PRIi64 ", %" PRIi64 " bp", nl, loci_bp);

    //////////////////////////////////////////////
    // (3) PAF pass 2, one query contig at a time (in byte order): departures
    //////////////////////////////////////////////

    int64_t nq = st->qc.n;
    int64_t *qorder = st_malloc(sizeof(int64_t) * (size_t)(nq > 0 ? nq : 1));
    {
        char **names = st_malloc(sizeof(char *) * (size_t)(nq > 0 ? nq : 1));
        for (int64_t i = 0; i < nq; i++) {
            qorder[i] = i;
            names[i] = st->qc.a[i].name;
        }
        qsort_names = names;
        qsort(qorder, (size_t)nq, sizeof(int64_t), cmp_qc_index);
        free(names);
    }
    // genomes of the query contigs
    for (int64_t i = 0; i < nq; i++) st->qc.a[i].genome = genome_index(st, st->qc.a[i].name);

    st->inside = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    st->contig_end_inside = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    st->maxstr = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    st->minstr = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    st->comp = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    st->mech = st_calloc(nl > 0 ? 3 * nl : 1, sizeof(int64_t));
    st->mask = st_calloc(nl > 0 ? nl : 1, sizeof(uint32_t));
    st->bubble = st_calloc(nl > 0 ? nl : 1, sizeof(bool));
    st->sv_detail = st_calloc(nl > 0 ? nl : 1, sizeof(CharVec));
    int64_t *rank_to_node = st_malloc(sizeof(int64_t) * (size_t)(nn > 0 ? nn : 1));
    for (int64_t i = 0; i < nn; i++) rank_to_node[N[i].rank] = i;
    int64_t *passno_stamp = st_malloc(sizeof(int64_t) * (size_t)(nl > 0 ? nl : 1));
    int64_t *passno_val = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    for (int64_t k = 0; k < nl; k++) passno_stamp[k] = -1;
    bool *multi_pass_locus = st_calloc(nl > 0 ? nl : 1, sizeof(bool));
    int64_t npass2 = 0;
    DRunVec druns = { 0 };
    for (int64_t k = 0; k < nq; k++) {
        load_contig_runs(st, in, &st->qc.a[qorder[k]], &sc, &druns);
        // passes record the contig by its rank, so that sorting passes by creation order sorts them by contig
        walk_contig(st, k, &druns, rank_to_node, passno_stamp, passno_val, multi_pass_locus, &npass2);
    }
    free(passno_stamp);
    free(passno_val);
    {
        CharVec m = { 0 };
        cv_printf(&m, "departures: ");
        bool first = 1;
        for (int k = 0; k < K_NUM; k++) {
            if (st->nkind[k]) {
                cv_printf(&m, "%s%s %" PRIi64, first ? "" : ", ", kind_names[k], st->nkind[k]);
                first = 0;
            }
        }
        ulog("%s", cv_str(&m));
        VEC_FREE(m);
    }

    //////////////////////////////////////////////
    // (4) Per-locus tallies
    //////////////////////////////////////////////

    // shared-alt: an alt node an inside excursion uses that is also used anywhere else
    sort_nt = &st->nt;
    if (st->altin.n > 0) qsort(st->altin.a, (size_t)st->altin.n, sizeof(AltIn), cmp_altin);
    {
        int64_t j = 0;
        for (int64_t k = 0; k < st->altin.n; k++) {
            if (j == 0 || st->altin.a[j - 1].node != st->altin.a[k].node || st->altin.a[j - 1].li != st->altin.a[k].li) {
                st->altin.a[j++] = st->altin.a[k];
            }
        }
        st->altin.n = j;
    }
    int64_t n_shared = 0, n_altin_nodes = 0, n_altin_out = 0, n_altin_multi = 0;
    for (int64_t k = 0; k < st->altin.n;) {
        int64_t j = k;
        while (j < st->altin.n && st->altin.a[j].node == st->altin.a[k].node) j++;
        Node *nd = &N[st->altin.a[k].node];
        n_altin_nodes++;
        if (nd->flags & NODE_ALT_OUT) n_altin_out++;
        if (j - k > 1) n_altin_multi++;
        if ((nd->flags & NODE_ALT_OUT) || j - k > 1) {
            nd->flags |= NODE_SHARED;
            n_shared++;
            for (int64_t g = k; g < j; g++) st->mask[st->altin.a[g].li] |= 1u << R_SHARED_ALT;
        }
        k = j;
    }
    // strings and compensating contigs
    for (int64_t li = 0; li < nl; li++) {
        st->maxstr[li] = st->minstr[li] = L[li].e - L[li].s + 2 * a->caf_trim;
    }
    if (st->passes.n > 0) qsort(st->passes.a, (size_t)st->passes.n, sizeof(Pass), cmp_pass);
    Pass *P = st->passes.a;
    int64_t np = st->passes.n;
    for (int64_t k = 0; k < np; k++) {
        int64_t sp = L[P[k].li].e - L[P[k].li].s;
        int64_t s_ = sp + P[k].ins - P[k].del + 2 * a->caf_trim;
        if (s_ > st->maxstr[P[k].li]) st->maxstr[P[k].li] = s_;
        if (s_ < st->minstr[P[k].li]) st->minstr[P[k].li] = s_;
        if (P[k].has_alt) st->bubble[P[k].li] = 1;
    }
    for (int64_t k = 0; k < np;) { // groups of (locus, contig)
        int64_t j = k;
        bool any_comp = 0, alt = 0, rep = 0;
        while (j < np && P[j].li == P[k].li && P[j].contig == P[k].contig) {
            if (P[j].ins > 0 && P[j].del > 0) {
                any_comp = 1;
                alt |= P[j].has_alt;
                rep |= P[j].has_replace;
            }
            j++;
        }
        if (any_comp) {
            st->comp[P[k].li]++;
            st->mech[3 * P[k].li + (alt ? 0 : rep ? 1 : 2)]++;
        }
        k = j;
    }
    {
        int64_t nloci2 = 0;
        for (int64_t li = 0; li < nl; li++) nloci2 += multi_pass_locus[li];
        ulog("contigs passing a locus more than once (each pass its own allele): %" PRIi64 " (locus, contig) pairs in %"
             PRIi64 " candidate loci", npass2, nloci2);
    }
    free(multi_pass_locus);
    for (int64_t li = 0; li < nl; li++) {
        if (st->maxstr[li] > a->max_string) st->mask[li] |= 1u << R_STRING;
    }

    //////////////////////////////////////////////
    // SV features of the bubble candidates' alt nodes
    //////////////////////////////////////////////

    uint32_t *sv_reasons = st_calloc(nl > 0 ? nl : 1, sizeof(uint32_t));
    VEC(SvRow) sv_rows = { 0 };
    bool sv_on = a->sv_veto || a->sv_report != NULL;
    if (sv_on) {
        double t_sv = now_seconds();
        int64_t nbloci = 0, nnodes = 0, nscan = 0;
        KmerVec tight = { 0 }, wide = { 0 }, qf = { 0 }, qr = { 0 };
        VEC(SvUse) uses = { 0 };
        I64Vec tmp = { 0 };
        for (int64_t k = 0; k < np;) {
            int64_t li = P[k].li, j = k;
            while (j < np && P[j].li == li) j++;
            if (!st->bubble[li]) {
                k = j;
                continue;
            }
            nbloci++;
            uses.n = 0;
            for (int64_t g = k; g < j; g++) {
                for (int64_t u = P[g].alt0; u < P[g].alt0 + P[g].nalt; u++) {
                    AltUse *au = &st->altuses.a[u];
                    SvUse x = { N[au->node].rank, au->node, g, au->ri, P[g].contig,
                                (au->rev ? -1 : 1) * P[g].dir > 0 };
                    VEC_PUSH(uses, x);
                }
            }
            if (uses.n > 0) qsort(uses.a, (size_t)uses.n, sizeof(SvUse), cmp_svuse);
            bool have_windows = 0;
            for (int64_t u = 0; u < uses.n;) {
                int64_t v = u, fwd = 0, rev = 0;
                tmp.n = 0;
                while (v < uses.n && uses.a[v].node == uses.a[u].node) {
                    if (v == u || uses.a[v].pass != uses.a[v - 1].pass || uses.a[v].ri != uses.a[v - 1].ri) {
                        if (uses.a[v].fwd) fwd++; else rev++;
                        VEC_PUSH(tmp, uses.a[v].contig);
                    }
                    v++;
                }
                Node *nd = &N[uses.a[u].node];
                u = v;
                int64_t len = (nd->flags & NODE_FASTA) ? nd->seq_len : 0;
                if (len < a->sv_min_bp) continue;
                nnodes++;
                // distinct contigs and genomes of its users
                qsort(tmp.a, (size_t)tmp.n, sizeof(int64_t), cmp_i64);
                int64_t ncont = 0;
                for (int64_t g = 0; g < tmp.n; g++) {
                    if (g == 0 || tmp.a[g] != tmp.a[g - 1]) tmp.a[ncont++] = tmp.a[g];
                }
                tmp.n = ncont;
                for (int64_t g = 0; g < tmp.n; g++) tmp.a[g] = st->qc.a[qorder[tmp.a[g]]].genome;
                qsort(tmp.a, (size_t)tmp.n, sizeof(int64_t), cmp_i64);
                int64_t ngen = 0;
                for (int64_t g = 0; g < tmp.n; g++) {
                    if (g == 0 || tmp.a[g] != tmp.a[g - 1]) ngen++;
                }
                char *hs = fwd >= rev ? stString_copy(nd->seq) : revcomp(nd->seq, len);
                char *hr = revcomp(hs, len);
                if (!have_windows) {
                    int64_t a0 = L[li].s - a->sv_tight > 0 ? L[li].s - a->sv_tight : 0;
                    int64_t a1 = L[li].e + a->sv_tight < st->G ? L[li].e + a->sv_tight : st->G;
                    kset(st->seq + a0, a1 > a0 ? a1 - a0 : 0, a->sv_k, &tight);
                    a0 = L[li].s - a->sv_wide > 0 ? L[li].s - a->sv_wide : 0;
                    a1 = L[li].e + a->sv_wide < st->G ? L[li].e + a->sv_wide : st->G;
                    kset(st->seq + a0, a1 > a0 ? a1 - a0 : 0, a->sv_k, &wide);
                    have_windows = 1;
                }
                kset(hs, len, a->sv_k, &qf);
                kset(hr, len, a->sv_k, &qr);
                double kf = contain(&qf, &tight), kfw = contain(&qf, &wide), krw = contain(&qr, &wide);
                bool novel = kf < a->sv_max_local;
                bool inv = krw >= a->sv_inv_min && krw >= kfw + a->sv_inv_margin;
                int64_t nontr = -1;
                if (a->sv_report != NULL || novel || inv) {
                    int64_t nc;
                    int64_t *c = tandem_scan(hs, len, a->min_period, a->max_period, a->window, a->min_ident, 1, &nc);
                    int64_t cov = 0;
                    for (int64_t g = 0; g < nc; g++) cov += c[3 * g + 1] - c[3 * g];
                    free(c);
                    nontr = len - cov;
                    nscan++;
                }
                uint32_t why = 0;
                if (fwd && rev) why |= 1u << R_SV_ORIENT;
                if (nontr >= 0 && nontr >= a->sv_min_bp) {
                    if (inv) why |= 1u << R_SV_INV;
                    if (novel) why |= 1u << R_SV_NOVEL;
                }
                const char *name = starts_with(nd->name, st->hub_t) ? nd->name + hub_t_len : nd->name;
                char nontr_s[32];
                if (nontr >= 0) snprintf(nontr_s, sizeof(nontr_s), "%" PRIi64, nontr);
                else snprintf(nontr_s, sizeof(nontr_s), ".");
                for (int r = R_SV_ORIENT; r <= R_SV_NOVEL; r++) {
                    if (!(why & (1u << r))) continue;
                    sv_reasons[li] |= 1u << r;
                    CharVec *d = &st->sv_detail[li];
                    cv_printf(d, "%s%s:%s:LN%" PRIi64 ":nonTR%s:kF%.2f:kFw%.2f:kRw%.2f:fwd%" PRIi64 "/rev%" PRIi64,
                              d->n ? " " : "", name, reason_names[r], len, nontr_s, kf, kfw, krw, fwd, rev);
                }
                if (a->sv_report != NULL) {
                    SvRow row = { li, len, fwd, rev, tmp.n, ngen, nontr, (char *)name, kf, kfw, krw, why };
                    VEC_PUSH(sv_rows, row);
                }
                free(hs);
                free(hr);
            }
            k = j;
        }
        if (a->sv_veto) {
            for (int64_t li = 0; li < nl; li++) st->mask[li] |= sv_reasons[li];
        }
        int64_t cnt[3] = { 0, 0, 0 };
        for (int64_t li = 0; li < nl; li++) {
            for (int r = 0; r < 3; r++) cnt[r] += (sv_reasons[li] >> (R_SV_ORIENT + r)) & 1;
        }
        ulog("SV features (%s): bubble candidate loci %" PRIi64 ", alt nodes >= %" PRIi64 " bp %" PRIi64
             ", tandem-scanned %" PRIi64 "; loci with an SV reason: sv-orient %" PRIi64 ", sv-inv %" PRIi64
             ", sv-novel %" PRIi64 "; %.1fs", a->sv_veto ? "veto on" : "veto OFF, report only", nbloci, a->sv_min_bp,
             nnodes, nscan, cnt[0], cnt[1], cnt[2], now_seconds() - t_sv);
        VEC_FREE(tight);
        VEC_FREE(wide);
        VEC_FREE(qf);
        VEC_FREE(qr);
        VEC_FREE(uses);
        VEC_FREE(tmp);
    }

    //////////////////////////////////////////////
    // Decide
    //////////////////////////////////////////////

    int64_t nkeep = 0;
    bool *kept = st_calloc(nl > 0 ? nl : 1, sizeof(bool));
    I64Vec keep = { 0 };
    int64_t first_count[R_NUM], gate_ok_count[R_NUM], gate_ok = 0;
    memset(first_count, 0, sizeof(first_count));
    memset(gate_ok_count, 0, sizeof(gate_ok_count));
    for (int64_t li = 0; li < nl; li++) {
        uint32_t m = st->mask[li];
        if (L[li].e - L[li].s > a->max_span) m |= 1u << R_SPAN;
        {
            int64_t i = lower_bound(st->uf_s.a, st->uf_s.n, L[li].e) - 1;
            if (i >= 0 && st->uf_emax.a[i] > L[li].s) m |= 1u << R_UNFRAMED;
        }
        bool gated = (a->gate == GATE_VARIABLE && st->inside[li] == 0) ||
                     (a->gate == GATE_COMP && st->comp[li] < a->min_haps);
        if (gated) m |= 1u << R_GATE;
        if (!a->sv_veto) m &= ~((1u << R_SV_ORIENT) | (1u << R_SV_INV) | (1u << R_SV_NOVEL));
        if (!a->keep_cigar_no_boundary && m == 0 && !st->bubble[li]) {
            int64_t j = lower_bound(st->bounds.a, st->bounds.n, L[li].s);
            bool has_boundary = j < st->bounds.n && st->bounds.a[j] <= L[li].e;
            if (!has_boundary) m = 1u << R_CIGAR_NOBOUNDARY;
        }
        st->mask[li] = m;
        if (m == 0) {
            kept[li] = 1;
            VEC_PUSH(keep, li);
            nkeep++;
        } else {
            first_count[first_reason(m)]++;
        }
        if (!(m & (1u << R_GATE))) {
            gate_ok++;
            for (int r = 0; r < R_NUM; r++) gate_ok_count[r] += (m >> r) & 1;
        }
    }
    {
        CharVec m = { 0 };
        cv_printf(&m, "gate %s minHaps %" PRIi64 ": %" PRIi64 " candidate loci pass the gate; vetoes among them (a locus "
                  "can have several): ", a->gate == GATE_COMP ? "comp" : a->gate == GATE_VARIABLE ? "variable" : "none",
                  a->min_haps, gate_ok);
        bool first = 1;
        for (int r = 0; r < R_NUM; r++) {
            if (r == R_GATE || !gate_ok_count[r]) continue;
            cv_printf(&m, "%s%s %" PRIi64, first ? "" : ", ", reason_names[r], gate_ok_count[r]);
            first = 0;
        }
        ulog("%s", cv_str(&m));
        m.n = 0;
        cv_printf(&m, "not cleared (first reason): ");
        first = 1;
        for (int r = 0; r < R_NUM; r++) {
            if (!first_count[r]) continue;
            cv_printf(&m, "%s%s %" PRIi64, first ? "" : ", ", reason_names[r], first_count[r]);
            first = 0;
        }
        ulog("%s", cv_str(&m));
        VEC_FREE(m);
    }
    if (n_shared) {
        ulog("shared-alt: %" PRIi64 " alt nodes used by an inside excursion and elsewhere", n_shared);
    }
    if (a->sv_veto) {
        int64_t nfirst = 0, nother = 0, nother_gate = 0;
        for (int64_t li = 0; li < nl; li++) {
            if (!st->mask[li]) continue;
            int fr = first_reason(st->mask[li]);
            if (fr >= R_SV_ORIENT && fr <= R_SV_NOVEL) nfirst++;
            else if (sv_reasons[li]) {
                nother++;
                if (fr == R_GATE) nother_gate++;
            }
        }
        ulog("svVeto: %" PRIi64 " loci not cleared because of an SV reason alone (first reason sv-*); %" PRIi64 " more "
             "loci with an SV reason were rejected for another reason first (gate %" PRIi64 ")", nfirst, nother,
             nother_gate);
        CharVec rs = { 0 };
        for (int64_t li = 0; li < nl; li++) {
            int fr = st->mask[li] ? first_reason(st->mask[li]) : -1;
            if (fr >= R_SV_ORIENT && fr <= R_SV_NOVEL) {
                const char *c;
                int64_t ls;
                to_local(st, L[li].s, &c, &ls);
                reasons_string(st->mask[li], &rs);
                ulog("svVeto: locus %" PRIi64 " %s:%" PRIi64 "-%" PRIi64 " period %" PRIi64 " compHaps %" PRIi64 " %s: %s",
                     li, c, ls, ls + L[li].e - L[li].s, L[li].p, st->comp[li], cv_str(&rs), cv_str(&st->sv_detail[li]));
            }
        }
        for (int64_t li = 0; li < nl; li++) {
            int fr = st->mask[li] ? first_reason(st->mask[li]) : -1;
            if (sv_reasons[li] && !(fr >= R_SV_ORIENT && fr <= R_SV_NOVEL) &&
                (sv_reasons[li] & ((1u << R_SV_ORIENT) | (1u << R_SV_INV)))) {
                const char *c;
                int64_t ls;
                to_local(st, L[li].s, &c, &ls);
                reasons_string(st->mask[li], &rs);
                ulog("svVeto (also rejected by %s): locus %" PRIi64 " %s:%" PRIi64 "-%" PRIi64 " %s: %s", reason_names[fr],
                     li, c, ls, ls + L[li].e - L[li].s, cv_str(&rs), cv_str(&st->sv_detail[li]));
            }
        }
        VEC_FREE(rs);
    }
    // alt nodes cleared whole: used only inside one kept locus
    int64_t n_clear_alt = 0;
    for (int64_t k = 0; k < st->altin.n;) {
        int64_t j = k;
        while (j < st->altin.n && st->altin.a[j].node == st->altin.a[k].node) j++;
        Node *nd = &N[st->altin.a[k].node];
        if (!(nd->flags & (NODE_ALT_OUT | NODE_SHARED)) && j - k == 1 && kept[st->altin.a[k].li]) {
            nd->flags |= NODE_CLEARED_ALT;
            nd->clear_locus = st->altin.a[k].li;
            n_clear_alt++;
        }
        k = j;
    }

    //////////////////////////////////////////////
    // Cleared hub intervals
    //////////////////////////////////////////////

    VEC(HubIv) hub = { 0 };
    {
        VEC(GNode) g = { 0 };
        for (int64_t i = 0; i < nn; i++) {
            if (N[i].base >= 0) {
                GNode x = { N[i].base, N[i].rank, i };
                VEC_PUSH(g, x);
            }
        }
        if (g.n > 0) qsort(g.a, (size_t)g.n, sizeof(GNode), cmp_gnode);
        I64Vec gs = { 0 };
        for (int64_t k = 0; k < g.n; k++) VEC_PUSH(gs, g.a[k].base);
        for (int64_t k = 0; k < keep.n; k++) {
            int64_t li = keep.a[k], s = L[li].s, e = L[li].e;
            int64_t j = upper_bound(gs.a, gs.n, s) - 1;
            if (j < 0) j = 0;
            while (j < g.n && g.a[j].base < e) {
                int64_t b = g.a[j].base, ln = N[g.a[j].node].ln;
                int64_t lo = s > b ? s : b, hi = e < b + ln ? e : b + ln;
                if (lo < hi) {
                    HubIv h = { li, g.a[j].node, lo - b, hi - b, 0 };
                    VEC_PUSH(hub, h);
                }
                j++;
            }
        }
        VEC_FREE(g);
        VEC_FREE(gs);
        for (int64_t i = 0; i < nn; i++) {
            if (N[i].flags & NODE_CLEARED_ALT) {
                HubIv h = { N[i].clear_locus, i, 0, N[i].tlen, 1 };
                VEC_PUSH(hub, h);
            }
        }
    }
    if (hub.n > 0) qsort(hub.a, (size_t)hub.n, sizeof(HubIv), cmp_hubiv_locus);
    int64_t *hub0 = st_calloc(nl + 1, sizeof(int64_t)); // hub intervals of locus li: [hub0[li], hub0[li+1])
    for (int64_t k = 0; k < hub.n; k++) hub0[hub.a[k].li + 1]++;
    for (int64_t li = 0; li < nl; li++) hub0[li + 1] += hub0[li];

    //////////////////////////////////////////////
    // (5) PAF pass 3: the contigs entering each kept locus
    //////////////////////////////////////////////

    I64Vec kstarts = { 0 };
    int64_t *keep_index = st_malloc(sizeof(int64_t) * (size_t)(nl > 0 ? nl : 1));
    for (int64_t k = 0; k < keep.n; k++) {
        VEC_PUSH(kstarts, L[keep.a[k]].s);
        keep_index[keep.a[k]] = k;
    }
    VEC(KC) enter = { 0 };
    if (keep.n > 0) {
        for (int64_t k = 0; k < nq; k++) {
            load_contig_runs(st, in, &st->qc.a[qorder[k]], &sc, &druns);
            for (int64_t j = 0; j < druns.n; j++) {
                DRun *r = &druns.a[j];
                if (r->o == 0) continue;
                for (int64_t g = upper_bound(kstarts.a, kstarts.n, r->rhi) - 1; g >= 0 && L[keep.a[g]].e > r->rlo; g--) {
                    KC x = { g, k };
                    VEC_PUSH(enter, x);
                }
            }
        }
        for (int64_t k = 0; k < np; k++) {
            if (kept[P[k].li]) {
                KC x = { keep_index[P[k].li], P[k].contig };
                VEC_PUSH(enter, x);
            }
        }
        if (enter.n > 0) qsort(enter.a, (size_t)enter.n, sizeof(KC), cmp_kc);
        int64_t j = 0;
        for (int64_t k = 0; k < enter.n; k++) {
            if (j == 0 || cmp_kc(&enter.a[j - 1], &enter.a[k]) != 0) enter.a[j++] = enter.a[k];
        }
        enter.n = j;
    }
    VEC_FREE(druns);

    //////////////////////////////////////////////
    // Write the BEDs
    //////////////////////////////////////////////

    // passes of kept loci by (locus, contig): P is sorted by (locus, contig, pass)
    int64_t *pass0 = st_calloc(nl + 1, sizeof(int64_t));
    for (int64_t k = 0; k < np; k++) pass0[P[k].li + 1]++;
    for (int64_t li = 0; li < nl; li++) pass0[li + 1] += pass0[li];

    int64_t nb_incl = 0, nb_strict = 0, kb = 0, nbub = 0, nbub_bp = 0, ncompalt = 0, ends_total = 0, ends_loci = 0;
    int64_t ends_noboundary = 0, hub_ref = 0, hub_alt = 0, max_str = 0;
    int64_t *nhaps_of = st_calloc(nl > 0 ? nl : 1, sizeof(int64_t));
    bool *boundary_of = st_calloc(nl > 0 ? nl : 1, sizeof(bool));
    FILE *loci_f = a->loci_bed ? open_out(a->loci_bed) : NULL;
    if (loci_f) fputs(loci_header, loci_f);
    {
        I64Vec strs = { 0 }, gens = { 0 };
        int64_t e0 = 0;
        for (int64_t k = 0; k < keep.n; k++) {
            int64_t li = keep.a[k], s = L[li].s, e = L[li].e, span = e - s;
            int64_t j = lower_bound(st->bounds.a, st->bounds.n, s);
            bool incl = j < st->bounds.n && st->bounds.a[j] <= e;
            int64_t j2 = upper_bound(st->bounds.a, st->bounds.n, s);
            bool strict = j2 < st->bounds.n && st->bounds.a[j2] < e;
            nb_incl += incl;
            nb_strict += strict;
            // the contigs entering the locus, their genomes and their pass strings
            while (e0 < enter.n && enter.a[e0].k < k) e0++;
            int64_t e1 = e0;
            while (e1 < enter.n && enter.a[e1].k == k) e1++;
            strs.n = gens.n = 0;
            int64_t pp = pass0[li];
            for (int64_t g = e0; g < e1; g++) {
                int64_t contig = enter.a[g].contig;
                VEC_PUSH(gens, st->qc.a[qorder[contig]].genome);
                while (pp < pass0[li + 1] && P[pp].contig < contig) pp++;
                if (pp < pass0[li + 1] && P[pp].contig == contig) {
                    while (pp < pass0[li + 1] && P[pp].contig == contig) {
                        VEC_PUSH(strs, span + P[pp].ins - P[pp].del);
                        pp++;
                    }
                } else {
                    VEC_PUSH(strs, span);
                }
            }
            e0 = e1;
            if (gens.n > 0) qsort(gens.a, (size_t)gens.n, sizeof(int64_t), cmp_i64);
            int64_t haps = 0;
            for (int64_t g = 0; g < gens.n; g++) {
                if (g == 0 || gens.a[g] != gens.a[g - 1]) haps++;
            }
            if (strs.n > 0) qsort(strs.a, (size_t)strs.n, sizeof(int64_t), cmp_i64);
            int64_t nalt = 0, nref = 0;
            for (int64_t h = hub0[li]; h < hub0[li + 1]; h++) {
                if (hub.a[h].alt) nalt++; else nref++;
            }
            hub_ref += nref;
            hub_alt += nalt;
            kb += span;
            if (st->bubble[li]) {
                nbub++;
                nbub_bp += span;
            }
            if (st->mech[3 * li] > 0) ncompalt++;
            ends_total += st->contig_end_inside[li];
            if (st->contig_end_inside[li]) {
                ends_loci++;
                if (!incl) ends_noboundary++;
            }
            if (st->maxstr[li] > max_str) max_str = st->maxstr[li];
            if (loci_f) {
                const char *c;
                int64_t ls;
                to_local(st, s, &c, &ls);
                fprintf(loci_f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%s\t%d\t%" PRIi64 "\t%" PRIi64
                        "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%"
                        PRIi64 "\t", c, ls, ls + span, li, L[li].p, st->bubble[li] ? "bubble" : "cigar", (int)incl, haps,
                        st->comp[li], st->inside[li], st->mech[3 * li], st->mech[3 * li + 1], st->mech[3 * li + 2],
                        st->minstr[li], st->maxstr[li], nref, nalt);
                for (int64_t g = 0; g < strs.n;) {
                    int64_t h = g;
                    while (h < strs.n && strs.a[h] == strs.a[g]) h++;
                    fprintf(loci_f, "%s%" PRIi64 "x%" PRIi64, g ? "," : "", strs.a[g], h - g);
                    g = h;
                }
                fprintf(loci_f, "\t%" PRIi64 "\n", st->contig_end_inside[li]);
            }
            // what the hub BED needs of each kept locus
            nhaps_of[li] = haps;
            boundary_of[li] = incl;
        }
        VEC_FREE(strs);
        VEC_FREE(gens);
    }
    if (loci_f) close_out(loci_f, a->loci_bed);
    if (a->hub_bed != NULL) {
        VEC(HubIv) hs = { 0 };
        VEC_GROW(hs, hub.n);
        if (hub.n > 0) memcpy(hs.a, hub.a, sizeof(HubIv) * (size_t)hub.n); // hub.a is NULL when nothing is cleared
        hs.n = hub.n;
        if (hs.n > 0) qsort(hs.a, (size_t)hs.n, sizeof(HubIv), cmp_hubiv_target);
        FILE *f = open_out(a->hub_bed);
        fputs(hub_header, f);
        for (int64_t k = 0; k < hs.n; k++) {
            HubIv *h = &hs.a[k];
            fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%s\t%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%"
                    PRIi64 "\n", N[h->node].name, h->ts, h->te, h->li, h->alt ? "alt" : "ref",
                    st->bubble[h->li] ? "bubble" : "cigar", (int64_t)boundary_of[h->li], nhaps_of[h->li],
                    st->comp[h->li], st->maxstr[h->li]);
        }
        close_out(f, a->hub_bed);
        VEC_FREE(hs);
    }
    if (a->veto_bed != NULL) {
        FILE *f = open_out(a->veto_bed);
        fputs(veto_header(a), f);
        CharVec rs = { 0 };
        for (int64_t li = 0; li < nl; li++) {
            if (kept[li]) continue;
            const char *c;
            int64_t ls;
            to_local(st, L[li].s, &c, &ls);
            reasons_string(st->mask[li], &rs);
            fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%s\t%s\t%" PRIi64, c, ls, ls + L[li].e - L[li].s,
                    L[li].p, reason_names[first_reason(st->mask[li])], cv_str(&rs), st->comp[li]);
            if (a->sv_veto) {
                fprintf(f, "\t%s", st->sv_detail[li].n ? cv_str(&st->sv_detail[li]) : "-");
            }
            fputc('\n', f);
        }
        VEC_FREE(rs);
        close_out(f, a->veto_bed);
    }
    if (a->sv_report != NULL) {
        FILE *f = open_out(a->sv_report);
        fputs(sv_header, f);
        CharVec rs = { 0 };
        for (int64_t k = 0; k < sv_rows.n; k++) {
            SvRow *r = &sv_rows.a[k];
            int64_t li = r->li;
            const char *c;
            int64_t ls;
            to_local(st, L[li].s, &c, &ls);
            rs.n = 0;
            cv_str(&rs);
            bool first = 1;
            for (int g = R_SV_ORIENT; g <= R_SV_NOVEL; g++) {
                if (r->why & (1u << g)) {
                    cv_printf(&rs, "%s%s", first ? "" : ",", reason_names[g]);
                    first = 0;
                }
            }
            char nontr_s[32];
            if (r->nontr >= 0) snprintf(nontr_s, sizeof(nontr_s), "%" PRIi64, r->nontr);
            else snprintf(nontr_s, sizeof(nontr_s), "None");
            fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%s\t%" PRIi64 "\t%s\t%" PRIi64 "\t%"
                    PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%" PRIi64 "\t%.4f\t%.4f\t%.4f\t%s\t%s\n", c, ls,
                    ls + L[li].e - L[li].s, li, L[li].p, kept[li] ? "cleared" : reason_names[first_reason(st->mask[li])],
                    st->comp[li], r->name, r->len, r->fwd, r->rev, r->ncontigs, r->ngenomes, r->kf, r->kfw, r->krw,
                    nontr_s, first ? "-" : cv_str(&rs));
        }
        VEC_FREE(rs);
        close_out(f, a->sv_report);
    }

    ulog("cleared loci %" PRIi64 ", %" PRIi64 " bp; class bubble %" PRIi64 " (%" PRIi64 " bp), cigar %" PRIi64
         "; loci with a comp haplotype through an alt node (critique altnode) %" PRIi64, nkeep, kb, nbub, nbub_bp,
         nkeep - nbub, ncompalt);
    ulog("boundary: rank-0 node boundary in [s,e] %" PRIi64 "; strictly inside (s<b<e, prototype) %" PRIi64, nb_incl,
         nb_strict);
    ulog("contig ends inside cleared loci: %" PRIi64 " in %" PRIi64 " loci (loci.bed contigEnds; such a locus also "
         "takes BAR's per-end path, boundary flag or not), of which without a rank-0 boundary %" PRIi64, ends_total,
         ends_loci, ends_noboundary);
    ulog("alt nodes cleared whole: %" PRIi64 " (alt nodes with inside excursions %" PRIi64 ", also used outside %"
         PRIi64 ", by several loci %" PRIi64 ")", n_clear_alt, n_altin_nodes, n_altin_out, n_altin_multi);
    ulog("hub intervals: %" PRIi64 " (%" PRIi64 " rank-0, %" PRIi64 " alt); max string between flanks %" PRIi64,
         hub.n, hub_ref, hub_alt, max_str);

    //////////////////////////////////////////////
    // (6) PAF pass 4, in input order: the cut
    //////////////////////////////////////////////

    // cleared intervals per target, merged (touching intervals too)
    Cleared *cl = st_calloc(nn > 0 ? nn : 1, sizeof(Cleared));
    {
        VEC(HubIv) hs = { 0 };
        VEC_GROW(hs, hub.n);
        if (hub.n > 0) memcpy(hs.a, hub.a, sizeof(HubIv) * (size_t)hub.n); // hub.a is NULL when nothing is cleared
        hs.n = hub.n;
        if (hs.n > 0) qsort(hs.a, (size_t)hs.n, sizeof(HubIv), cmp_hubiv_target);
        for (int64_t k = 0; k < hs.n;) {
            int64_t j = k;
            while (j < hs.n && hs.a[j].node == hs.a[k].node) j++;
            Cleared *c = &cl[hs.a[k].node];
            c->s = st_malloc(sizeof(int64_t) * (size_t)(j - k));
            c->e = st_malloc(sizeof(int64_t) * (size_t)(j - k));
            for (int64_t g = k; g < j; g++) {
                if (hs.a[g].te <= hs.a[g].ts) continue;
                if (c->n > 0 && hs.a[g].ts <= c->e[c->n - 1]) {
                    if (hs.a[g].te > c->e[c->n - 1]) c->e[c->n - 1] = hs.a[g].te;
                } else {
                    c->s[c->n] = hs.a[g].ts;
                    c->e[c->n] = hs.a[g].te;
                    c->n++;
                }
            }
            k = j;
        }
        VEC_FREE(hs);
    }
    // query sequences: R (matched intervals of changed records) and K (surviving matched intervals)
    stHash *rq_index = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
    VEC(char *) rq_names = { 0 };
    VEC(IvVec) rq_rem = { 0 }, rq_kept = { 0 };
    uint8_t *changed = st_calloc((size_t)(nrec / 8 + 1), 1);
    int64_t totc = 0, remr = 0, rema = 0, totnr = 0, remnr = 0;
    int64_t c_in = 0, c_unchanged = 0, c_split = 0, c_frags = 0, c_dropped = 0;
    FILE *out = a->out_paf == NULL ? stdout : open_out(a->out_paf);
    {
        if (fseeko(in, 0, SEEK_SET) != 0) st_errnoAbort("Could not rewind %s", a->paf);
        char *buf = NULL;
        size_t cap = 0;
        ssize_t n;
        int64_t ri = 0;
        FragVec frags = { 0 };
        IvVec removed = { 0 }, keptv = { 0 };
        while ((n = getline(&buf, &cap, in)) >= 0) {
            if (is_blank_line(buf, n)) continue;
            c_in++;
            int64_t len = n;
            if (buf[len - 1] == '\n') len--;
            VEC_GROW(sc.line, len + 1);
            memcpy(sc.line.a, buf, (size_t)len);
            sc.line.a[len] = '\0';
            split_fields(sc.line.a, &sc.fl);
            Rec r;
            parse_record(&sc.fl, &r, &sc.ops);
            int64_t node = node_get(&st->nt, r.t, 0);
            bool isref = starts_with(r.q, st->ref_q);
            Cleared *c = &cl[node];
            // matched columns, and those on cleared positions (log)
            record_runs(&r, &sc.ops, &sc.runs);
            for (int64_t k = 0; k < sc.runs.n; k++) {
                int64_t t0 = sc.runs.a[k].t0, rn = sc.runs.a[k].n, rm = 0;
                totc += rn;
                if (!isref) totnr += rn;
                if (c->n == 0) continue;
                for (int64_t j = upper_bound(c->s, c->n, t0 + rn) - 1; j >= 0 && c->e[j] > t0; j--) {
                    int64_t lo = t0 > c->s[j] ? t0 : c->s[j], hi = t0 + rn < c->e[j] ? t0 + rn : c->e[j];
                    if (hi > lo) rm += hi - lo;
                }
                if (N[node].flags & NODE_CLEARED_ALT) rema += rm; else remr += rm;
                if (!isref) remnr += rm;
            }
            bool cut = 0;
            if (!a->dry_run && c->n > 0) {
                removed.n = keptv.n = 0;
                cut = cut_record(&r, &sc.ops, c, &frags, &removed, &keptv);
            }
            if (!cut) {
                fwrite(buf, 1, (size_t)n, out);
                if (buf[n - 1] != '\n') fputc('\n', out);
                c_unchanged++;
                ri++;
                continue;
            }
            changed[ri / 8] |= (uint8_t)(1u << (ri % 8));
            ri++;
            if (frags.n == 0) {
                c_dropped++;
            } else {
                c_split++;
                c_frags += frags.n;
            }
            if (frags.n > 0) {
                int64_t n_eq = 0, n_m = 0;
                for (int64_t k = 0; k < sc.ops.n; k++) {
                    if (sc.ops.a[k].op == '=') n_eq += sc.ops.a[k].len;
                    if (sc.ops.a[k].op == 'M') n_m += sc.ops.a[k].len;
                }
                double m_ratio = n_m ? (double)(r.nmatch - n_eq) / (double)n_m : 0.0;
                if (m_ratio < 0.0) m_ratio = 0.0;
                if (m_ratio > 1.0) m_ratio = 1.0;
                for (int64_t k = 0; k < frags.n; k++) {
                    write_fragment(out, &sc.fl, &r, &frags.a[k], m_ratio);
                }
            }
            if (a->query_bed != NULL) {
                void *x = stHash_search(rq_index, r.q);
                int64_t qi;
                if (x == NULL) {
                    char *nm = stString_copy(r.q);
                    VEC_PUSH(rq_names, nm);
                    IvVec e1 = { 0 }, e2 = { 0 };
                    VEC_PUSH(rq_rem, e1);
                    VEC_PUSH(rq_kept, e2);
                    qi = rq_names.n - 1;
                    stHash_insert(rq_index, stString_copy(nm), (void *)(intptr_t)(qi + 1));
                } else {
                    qi = (int64_t)(intptr_t)x - 1;
                }
                IvVec *rv = &rq_rem.a[qi], *kv = &rq_kept.a[qi];
                VEC_GROW(*rv, removed.n);
                if (removed.n > 0) memcpy(rv->a + rv->n, removed.a, sizeof(Iv) * (size_t)removed.n);
                rv->n += removed.n;
                VEC_GROW(*kv, keptv.n);
                if (keptv.n > 0) memcpy(kv->a + kv->n, keptv.a, sizeof(Iv) * (size_t)keptv.n); // a dropped record keeps none
                kv->n += keptv.n;
            }
        }
        if (ferror(in)) st_errnoAbort("Error reading %s", a->paf);
        free(buf);
        for (int64_t k = 0; k < frags.cap; k++) VEC_FREE(frags.a[k].ops);
        free(frags.a);
        VEC_FREE(removed);
        VEC_FREE(keptv);
    }
    if (a->out_paf != NULL) close_out(out, a->out_paf);
    ulog("matched columns: all genomes %" PRIi64 ", removed on rank-0 nodes %" PRIi64 ", on alt nodes %" PRIi64
         " (%.2f%%); non-reference genomes %" PRIi64 ", removed %" PRIi64 " (%.2f%%)", totc, remr, rema,
         100.0 * (double)(remr + rema) / (double)(totc > 1 ? totc : 1), totnr, remnr,
         100.0 * (double)remnr / (double)(totnr > 1 ? totnr : 1));
    ulog("time %.1fs", now_seconds() - T0);
    ulog("cut_paf: records in %" PRIi64 ", unchanged %" PRIi64 ", cut %" PRIi64 " into %" PRIi64 " fragments, dropped %"
         PRIi64 ", records out %" PRIi64, c_in, c_unchanged, c_split, c_frags, c_dropped, c_unchanged + c_frags);

    //////////////////////////////////////////////
    // (7) PAF pass 5: the query BED (bases with a matched column in the input and none in the output)
    //////////////////////////////////////////////

    if (a->query_bed != NULL) {
        int64_t nqr = rq_names.n;
        for (int64_t i = 0; i < nqr; i++) merge_ivs(&rq_rem.a[i]);
        if (nqr > 0) {
            // unchanged records of these queries keep all their matched columns
            if (fseeko(in, 0, SEEK_SET) != 0) st_errnoAbort("Could not rewind %s", a->paf);
            char *buf = NULL;
            size_t cap = 0;
            ssize_t n;
            int64_t ri = 0;
            while ((n = getline(&buf, &cap, in)) >= 0) {
                if (is_blank_line(buf, n)) continue;
                int64_t my = ri++;
                if (changed[my / 8] & (1u << (my % 8))) continue;
                if (buf[n - 1] == '\n') buf[n - 1] = '\0';
                char *tab = strchr(buf, '\t');
                if (tab == NULL) continue;
                *tab = '\0';
                void *x = stHash_search(rq_index, buf);
                if (x == NULL) continue;
                *tab = '\t';
                IvVec *rv = &rq_rem.a[(int64_t)(intptr_t)x - 1], *kv = &rq_kept.a[(int64_t)(intptr_t)x - 1];
                split_fields(buf, &sc.fl);
                Rec r;
                parse_record(&sc.fl, &r, &sc.ops);
                int64_t q = 0;
                for (int64_t k = 0; k < sc.ops.n; k++) {
                    char op = sc.ops.a[k].op;
                    if (op == 'M' || op == '=' || op == 'X') {
                        Iv v;
                        int64_t ln = sc.ops.a[k].len;
                        if (r.rev) { v.s = r.qe - q - ln; v.e = r.qe - q; } else { v.s = r.qs + q; v.e = r.qs + q + ln; }
                        // only what can intersect R matters
                        int64_t lo = 0, hi = rv->n;
                        while (lo < hi) { // first R interval with e > v.s
                            int64_t mid = lo + (hi - lo) / 2;
                            if (rv->a[mid].e <= v.s) lo = mid + 1; else hi = mid;
                        }
                        if (lo < rv->n && rv->a[lo].s < v.e) {
                            // touching the last one (the record's previous op): extend it, as add_query_iv does
                            Iv *l = kv->n > 0 ? &kv->a[kv->n - 1] : NULL;
                            if (l != NULL && v.s == l->e) l->e = v.e;
                            else if (l != NULL && v.e == l->s) l->s = v.s;
                            else VEC_PUSH(*kv, v);
                        }
                        q += ln;
                    } else if (op == 'I') {
                        q += sc.ops.a[k].len;
                    }
                }
            }
            free(buf);
            if (ferror(in)) st_errnoAbort("Error reading %s", a->paf);
        }
        int64_t *qo = st_malloc(sizeof(int64_t) * (size_t)(nqr > 0 ? nqr : 1));
        for (int64_t i = 0; i < nqr; i++) qo[i] = i;
        qsort_names = rq_names.a;
        qsort(qo, (size_t)nqr, sizeof(int64_t), cmp_qc_index);
        FILE *f = open_out(a->query_bed);
        int64_t nqi = 0, nqb = 0;
        for (int64_t k = 0; k < nqr; k++) {
            int64_t i = qo[k];
            IvVec *rv = &rq_rem.a[i], *kv = &rq_kept.a[i];
            merge_ivs(kv);
            int64_t j = 0;
            for (int64_t g = 0; g < rv->n; g++) {
                int64_t cur = rv->a[g].s, e = rv->a[g].e;
                while (j < kv->n && kv->a[j].e <= cur) j++;
                for (int64_t h = j; h < kv->n && kv->a[h].s < e; h++) {
                    if (kv->a[h].s > cur) {
                        fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\n", rq_names.a[i], cur, kv->a[h].s);
                        nqi++;
                        nqb += kv->a[h].s - cur;
                    }
                    if (kv->a[h].e > cur) cur = kv->a[h].e;
                    if (cur >= e) break;
                }
                if (cur < e) {
                    fprintf(f, "%s\t%" PRIi64 "\t%" PRIi64 "\n", rq_names.a[i], cur, e);
                    nqi++;
                    nqb += e - cur;
                }
            }
            VEC_FREE(*rv);
            VEC_FREE(*kv);
        }
        close_out(f, a->query_bed);
        free(qo);
        ulog("cut_paf: query bases that lost every hub anchor: %" PRIi64 " bp in %" PRIi64 " intervals over %" PRIi64
             " query sequences", nqb, nqi, nqr);
    }
    fclose(in);
    write_log_file(a);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    for (int64_t i = 0; i < rq_names.n; i++) free(rq_names.a[i]);
    VEC_FREE(rq_names);
    VEC_FREE(rq_rem);
    VEC_FREE(rq_kept);
    stHash_destruct(rq_index);
    free(changed);
    for (int64_t i = 0; i < nn; i++) {
        free(cl[i].s);
        free(cl[i].e);
    }
    free(cl);
    free(pass0);
    free(nhaps_of);
    free(boundary_of);
    VEC_FREE(enter);
    VEC_FREE(kstarts);
    free(keep_index);
    free(hub0);
    VEC_FREE(hub);
    VEC_FREE(keep);
    free(kept);
    VEC_FREE(sv_rows);
    free(sv_reasons);
    for (int64_t li = 0; li < nl; li++) VEC_FREE(st->sv_detail[li]);
    free(st->sv_detail);
    free(rank_to_node);
    free(qorder);
    free(st->inside);
    free(st->contig_end_inside);
    free(st->maxstr);
    free(st->minstr);
    free(st->comp);
    free(st->mech);
    free(st->mask);
    free(st->bubble);
    VEC_FREE(st->passes);
    VEC_FREE(st->altuses);
    VEC_FREE(st->altin);
    VEC_FREE(st->loci);
    VEC_FREE(st->lstarts);
    VEC_FREE(st->bounds);
    VEC_FREE(st->uf_s);
    VEC_FREE(st->uf_emax);
    free(st->seq);
    free(st->contigs);
    free(st->clen);
    free(st->coff);
    free(corder);
    free(cpos);
    free_inputs(st, &sc, &ref_names, &ref_lens, ref_len_index);
    st_logInfo("Paffy unanchor is done, %.1f seconds have elapsed\n", now_seconds() - T0);
    return 0;
}
