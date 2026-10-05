/*
 * Unit tests for individual paf.h API functions.
 */

#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"

/* ---- helpers ---- */

/* paf_parse modifies its input (strtok_r), so always pass a copy */
static Paf *parse_str(const char *s, bool cigar) {
    char *copy = stString_copy(s);
    Paf *p = paf_parse(copy, cigar);
    free(copy);
    return p;
}

/* Build a PAF record programmatically */
static Paf *make_paf(const char *qname, int64_t qlen, int64_t qs, int64_t qe,
                     bool same_strand,
                     const char *tname, int64_t tlen, int64_t ts, int64_t te,
                     int64_t nm, int64_t nb, int64_t mq,
                     const char *cigar_str) {
    Paf *p = st_calloc(1, sizeof(Paf));
    p->query_name   = stString_copy(qname);
    p->query_length = qlen;
    p->query_start  = qs;
    p->query_end    = qe;
    p->target_name   = stString_copy(tname);
    p->target_length = tlen;
    p->target_start  = ts;
    p->target_end    = te;
    p->same_strand   = same_strand;
    p->num_matches   = nm;
    p->num_bases     = nb;
    p->mapping_quality = mq;
    p->tile_level = -1;
    p->chain_id   = -1;
    p->chain_score = -1;
    if (cigar_str) {
        char *cs = stString_copy(cigar_str);
        p->cigar = cigar_parse(cs);
        free(cs);
    }
    return p;
}

/* ---- 1. Cigar parsing ---- */

static void test_cigar_parse_empty(CuTest *tc) {
    Cigar *c = cigar_parse("");
    CuAssertTrue(tc, c == NULL);
}

static void test_cigar_parse_single(CuTest *tc) {
    char s[] = "10M";
    Cigar *c = cigar_parse(s);
    CuAssertTrue(tc, c != NULL);
    CuAssertIntEquals(tc, 1, cigar_count(c));
    CuAssertIntEquals(tc, match, cigar_get(c, 0)->op);
    CuAssertTrue(tc, cigar_get(c, 0)->length == 10);
    cigar_destruct(c);
}

static void test_cigar_parse_all_ops(CuTest *tc) {
    char s[] = "5M3I2D4=1X";
    Cigar *c = cigar_parse(s);
    CuAssertTrue(tc, c != NULL);
    CuAssertIntEquals(tc, 5, cigar_count(c));
    CuAssertIntEquals(tc, match,             cigar_get(c, 0)->op);
    CuAssertTrue(tc, cigar_get(c, 0)->length == 5);
    CuAssertIntEquals(tc, query_insert,      cigar_get(c, 1)->op);
    CuAssertTrue(tc, cigar_get(c, 1)->length == 3);
    CuAssertIntEquals(tc, query_delete,      cigar_get(c, 2)->op);
    CuAssertTrue(tc, cigar_get(c, 2)->length == 2);
    CuAssertIntEquals(tc, sequence_match,    cigar_get(c, 3)->op);
    CuAssertTrue(tc, cigar_get(c, 3)->length == 4);
    CuAssertIntEquals(tc, sequence_mismatch, cigar_get(c, 4)->op);
    CuAssertTrue(tc, cigar_get(c, 4)->length == 1);
    cigar_destruct(c);
}

static void test_cigar_parse_large_length(CuTest *tc) {
    char s[] = "1000000M";
    Cigar *c = cigar_parse(s);
    CuAssertTrue(tc, c != NULL);
    CuAssertIntEquals(tc, 1, cigar_count(c));
    CuAssertTrue(tc, cigar_get(c, 0)->length == 1000000);
    cigar_destruct(c);
}

/* ---- 2. Cigar accessors ---- */

static void test_cigar_count_get(CuTest *tc) {
    char s[] = "3M2I";
    Cigar *c = cigar_parse(s);
    CuAssertIntEquals(tc, 2, cigar_count(c));
    CuAssertIntEquals(tc, match,        cigar_get(c, 0)->op);
    CuAssertTrue(tc, cigar_get(c, 0)->length == 3);
    CuAssertIntEquals(tc, query_insert, cigar_get(c, 1)->op);
    CuAssertTrue(tc, cigar_get(c, 1)->length == 2);
    /* NULL -> 0 */
    CuAssertIntEquals(tc, 0, cigar_count(NULL));
    cigar_destruct(c);
}

/* ---- 3. PAF parsing ---- */

static void test_paf_parse_minimal(CuTest *tc) {
    Paf *paf = parse_str(
        "query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t50\t50\t255",
        true);
    CuAssertStrEquals(tc, "query1",  paf->query_name);
    CuAssertTrue(tc, paf->query_length == 100);
    CuAssertTrue(tc, paf->query_start  == 0);
    CuAssertTrue(tc, paf->query_end    == 50);
    CuAssertStrEquals(tc, "target1", paf->target_name);
    CuAssertTrue(tc, paf->target_length == 200);
    CuAssertTrue(tc, paf->target_start  == 10);
    CuAssertTrue(tc, paf->target_end    == 60);
    CuAssertTrue(tc, paf->num_matches   == 50);
    CuAssertTrue(tc, paf->num_bases     == 50);
    CuAssertTrue(tc, paf->mapping_quality == 255);
    CuAssertTrue(tc, paf->same_strand == true);
    CuAssertTrue(tc, paf->cigar        == NULL);
    CuAssertTrue(tc, paf->cigar_string == NULL);
    paf_destruct(paf);
}

static void test_paf_parse_with_cigar(CuTest *tc) {
    /* query span: 5+3=8, target span: 5+2=7 */
    Paf *paf = parse_str(
        "q1\t100\t0\t8\t+\tt1\t200\t0\t7\t8\t10\t60\tcg:Z:5M3I2D",
        true);
    CuAssertTrue(tc, paf->cigar        != NULL);
    CuAssertTrue(tc, paf->cigar_string == NULL);
    CuAssertIntEquals(tc, 3, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, match,        cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    CuAssertIntEquals(tc, query_insert, cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 3);
    CuAssertIntEquals(tc, query_delete, cigar_get(paf->cigar, 2)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 2)->length == 2);
    paf_destruct(paf);
}

static void test_paf_parse_cigar_string_mode(CuTest *tc) {
    Paf *paf = parse_str(
        "q1\t100\t0\t8\t+\tt1\t200\t0\t7\t8\t10\t60\tcg:Z:5M3I2D",
        false);
    CuAssertTrue(tc, paf->cigar        == NULL);
    CuAssertTrue(tc, paf->cigar_string != NULL);
    CuAssertStrEquals(tc, "5M3I2D", paf->cigar_string);
    paf_destruct(paf);
}

static void test_paf_parse_optional_tags(CuTest *tc) {
    Paf *paf = parse_str(
        "q1\t100\t0\t50\t+\tt1\t200\t0\t50\t50\t50\t60\t"
        "tp:A:P\tAS:i:42\ttl:i:2\tcn:i:5\ts1:i:100",
        true);
    CuAssertTrue(tc, paf->type        == 'P');
    CuAssertTrue(tc, paf->score       == 42);
    CuAssertTrue(tc, paf->tile_level  == 2);
    CuAssertTrue(tc, paf->chain_id    == 5);
    CuAssertTrue(tc, paf->chain_score == 100);
    /* absent optional tags default to -1 */
    paf_destruct(paf);
}

static void test_paf_parse_strand(CuTest *tc) {
    Paf *pos = parse_str(
        "q1\t100\t0\t50\t+\tt1\t200\t0\t50\t50\t50\t60", true);
    CuAssertTrue(tc, pos->same_strand == true);
    paf_destruct(pos);

    Paf *neg = parse_str(
        "q1\t100\t0\t50\t-\tt1\t200\t0\t50\t50\t50\t60", true);
    CuAssertTrue(tc, neg->same_strand == false);
    paf_destruct(neg);
}

/* ---- 4. Roundtrip ---- */

static void test_paf_roundtrip_no_cigar(CuTest *tc) {
    /* Parse, print, re-parse, re-print: second and third strings must match */
    Paf *paf1 = parse_str(
        "query1\t100\t0\t50\t+\ttarget1\t200\t10\t60\t50\t50\t255",
        true);
    char *s1      = paf_print(paf1);
    char *s1_copy = stString_copy(s1);
    free(s1);

    Paf *paf2 = parse_str(s1_copy, true);
    char *s2  = paf_print(paf2);

    CuAssertStrEquals(tc, s1_copy, s2);
    free(s1_copy);
    free(s2);
    paf_destruct(paf1);
    paf_destruct(paf2);
}

static void test_paf_roundtrip_with_cigar(CuTest *tc) {
    /* 5M3I2D: query span=8, target span=7 */
    Paf *paf1 = parse_str(
        "q1\t100\t0\t8\t+\tt1\t200\t0\t7\t8\t10\t60\tcg:Z:5M3I2D",
        true);
    char *s1      = paf_print(paf1);
    char *s1_copy = stString_copy(s1);
    free(s1);

    Paf *paf2 = parse_str(s1_copy, true);
    char *s2  = paf_print(paf2);

    CuAssertStrEquals(tc, s1_copy, s2);
    CuAssertIntEquals(tc, 3, cigar_count(paf2->cigar));
    CuAssertIntEquals(tc, match,        cigar_get(paf2->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf2->cigar, 0)->length == 5);
    CuAssertIntEquals(tc, query_insert, cigar_get(paf2->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf2->cigar, 1)->length == 3);
    CuAssertIntEquals(tc, query_delete, cigar_get(paf2->cigar, 2)->op);
    CuAssertTrue(tc, cigar_get(paf2->cigar, 2)->length == 2);

    free(s1_copy);
    free(s2);
    paf_destruct(paf1);
    paf_destruct(paf2);
}

/* ---- 5. File I/O ---- */

static void test_paf_read_write(CuTest *tc) {
    FILE *fh = tmpfile();
    CuAssertTrue(tc, fh != NULL);

    fprintf(fh, "q1\t100\t0\t50\t+\tt1\t200\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q2\t200\t10\t60\t-\tt2\t300\t20\t70\t50\t50\t30\n");
    fprintf(fh, "q3\t150\t5\t55\t+\tt3\t250\t15\t65\t50\t50\t40\n");
    rewind(fh);

    Paf *p1   = paf_read2(fh);
    Paf *p2   = paf_read2(fh);
    Paf *p3   = paf_read2(fh);
    Paf *pend = paf_read2(fh);

    CuAssertTrue(tc, p1   != NULL);
    CuAssertTrue(tc, p2   != NULL);
    CuAssertTrue(tc, p3   != NULL);
    CuAssertTrue(tc, pend == NULL);

    CuAssertStrEquals(tc, "q1", p1->query_name);
    CuAssertTrue(tc, p1->query_length == 100);
    CuAssertTrue(tc, p1->same_strand  == true);

    CuAssertStrEquals(tc, "q2", p2->query_name);
    CuAssertTrue(tc, p2->same_strand == false);

    CuAssertStrEquals(tc, "q3", p3->query_name);
    CuAssertTrue(tc, p3->query_start == 5);

    paf_destruct(p1);
    paf_destruct(p2);
    paf_destruct(p3);
    fclose(fh);
}

static void test_read_write_pafs_list(CuTest *tc) {
    stList *out = stList_construct3(0, (void(*)(void*))paf_destruct);
    stList_append(out, make_paf("qa", 100, 0, 50, true,  "ta", 200, 0, 50, 50, 50, 60, NULL));
    stList_append(out, make_paf("qb", 100, 0, 50, false, "tb", 200, 0, 50, 50, 50, 60, NULL));
    stList_append(out, make_paf("qc", 100, 0, 50, true,  "tc", 200, 0, 50, 50, 50, 60, NULL));

    FILE *fh = tmpfile();
    CuAssertTrue(tc, fh != NULL);
    write_pafs(fh, out);
    rewind(fh);
    stList *in = read_pafs(fh, false);
    fclose(fh);

    CuAssertIntEquals(tc, 3, stList_length(in));
    CuAssertStrEquals(tc, "qa", ((Paf*)stList_get(in, 0))->query_name);
    CuAssertTrue(tc, ((Paf*)stList_get(in, 0))->same_strand == true);
    CuAssertStrEquals(tc, "qb", ((Paf*)stList_get(in, 1))->query_name);
    CuAssertTrue(tc, ((Paf*)stList_get(in, 1))->same_strand == false);
    CuAssertStrEquals(tc, "qc", ((Paf*)stList_get(in, 2))->query_name);

    stList_destruct(out);
    stList_destruct(in);
}

/* ---- 6. PAF Stats ---- */

static void test_paf_stats_calc_all_match(CuTest *tc) {
    Paf *paf = make_paf("q", 100, 0, 10, true, "t", 100, 0, 10, 10, 10, 60, "10M");
    int64_t mat=0, mis=0, qi=0, qd=0, qib=0, qdb=0;
    paf_stats_calc(paf, &mat, &mis, &qi, &qd, &qib, &qdb, true);
    CuAssertTrue(tc, mat == 10);
    CuAssertTrue(tc, mis == 0 && qi == 0 && qd == 0 && qib == 0 && qdb == 0);
    paf_destruct(paf);
}

static void test_paf_stats_calc_mixed(CuTest *tc) {
    /* "3=2X1I2D": query span=3+2+1=6, target span=3+2+2=7 */
    Paf *paf = make_paf("q", 100, 0, 6, true, "t", 100, 0, 7, 5, 8, 60, "3=2X1I2D");
    int64_t mat=0, mis=0, qi=0, qd=0, qib=0, qdb=0;
    paf_stats_calc(paf, &mat, &mis, &qi, &qd, &qib, &qdb, true);
    CuAssertTrue(tc, mat == 3);
    CuAssertTrue(tc, mis == 2);
    CuAssertTrue(tc, qi  == 1 && qib == 1);
    CuAssertTrue(tc, qd  == 1 && qdb == 2);
    paf_destruct(paf);
}

static void test_paf_stats_calc_zero_flag(CuTest *tc) {
    Paf *paf = make_paf("q", 100, 0, 5, true, "t", 100, 0, 5, 5, 5, 60, "5M");
    int64_t mat=0, mis=0, qi=0, qd=0, qib=0, qdb=0;

    /* accumulate twice without zeroing */
    paf_stats_calc(paf, &mat, &mis, &qi, &qd, &qib, &qdb, false);
    paf_stats_calc(paf, &mat, &mis, &qi, &qd, &qib, &qdb, false);
    CuAssertTrue(tc, mat == 10);

    /* zero_counts=true resets before accumulating */
    paf_stats_calc(paf, &mat, &mis, &qi, &qd, &qib, &qdb, true);
    CuAssertTrue(tc, mat == 5);

    paf_destruct(paf);
}

/* ---- 7. PAF Invert ---- */

static void test_paf_invert_same_strand(CuTest *tc) {
    /* 5M3I2D: query span=5+3=8, target span=5+2=7 */
    Paf *paf = make_paf("query", 100, 10, 18, true, "target", 200, 20, 27, 8, 10, 60, "5M3I2D");
    paf_invert(paf);

    /* names swapped */
    CuAssertStrEquals(tc, "target", paf->query_name);
    CuAssertStrEquals(tc, "query",  paf->target_name);
    /* coords swapped */
    CuAssertTrue(tc, paf->query_start  == 20);
    CuAssertTrue(tc, paf->query_end    == 27);
    CuAssertTrue(tc, paf->query_length == 200);
    CuAssertTrue(tc, paf->target_start  == 10);
    CuAssertTrue(tc, paf->target_end    == 18);
    CuAssertTrue(tc, paf->target_length == 100);
    /* same_strand unchanged */
    CuAssertTrue(tc, paf->same_strand == true);
    /* I<->D swapped, order unchanged (same_strand) → 5M3D2I */
    CuAssertIntEquals(tc, 3, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, match,        cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    CuAssertIntEquals(tc, query_delete, cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 3);
    CuAssertIntEquals(tc, query_insert, cigar_get(paf->cigar, 2)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 2)->length == 2);

    paf_destruct(paf);
}

static void test_paf_invert_opposite_strand(CuTest *tc) {
    /* 5M3I: query span=5+3=8, target span=5 */
    Paf *paf = make_paf("query", 100, 10, 18, false, "target", 200, 20, 25, 5, 8, 60, "5M3I");
    paf_invert(paf);

    CuAssertTrue(tc, paf->same_strand == false);
    CuAssertStrEquals(tc, "target", paf->query_name);
    CuAssertStrEquals(tc, "query",  paf->target_name);

    /* I<->D swapped and then reversed (opposite strand): 5M3D → reversed → 3D5M */
    CuAssertIntEquals(tc, 2, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, query_delete, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 3);
    CuAssertIntEquals(tc, match,        cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 5);

    paf_destruct(paf);
}

static void test_paf_invert_double(CuTest *tc) {
    /* Double-invert must return to original state */
    Paf *paf = make_paf("query", 100, 10, 18, true, "target", 200, 20, 27, 8, 10, 60, "5M3I2D");
    char *orig = paf_print(paf);
    paf_invert(paf);
    paf_invert(paf);
    char *trip = paf_print(paf);
    CuAssertStrEquals(tc, orig, trip);
    free(orig);
    free(trip);
    paf_destruct(paf);
}

/* ---- 8. Aligned base count ---- */

static void test_aligned_bases(CuTest *tc) {
    /* "5M3I2D4=1X": M+=X count, I/D excluded → 5+4+1=10 */
    /* query span=5+3+4+1=13, target span=5+2+4+1=12 */
    Paf *paf = make_paf("q", 100, 0, 13, true, "t", 100, 0, 12, 10, 15, 60, "5M3I2D4=1X");
    CuAssertTrue(tc, paf_get_number_of_aligned_bases(paf) == 10);
    paf_destruct(paf);
}

/* ---- 9. Trimming ---- */

static void test_paf_trim_ends_zero(CuTest *tc) {
    Paf *paf = make_paf("q", 100, 5, 15, true, "t", 100, 5, 15, 10, 10, 60, "10M");
    paf_trim_ends(paf, 0);
    CuAssertTrue(tc, paf->query_start  == 5);
    CuAssertTrue(tc, paf->query_end    == 15);
    CuAssertTrue(tc, paf->target_start == 5);
    CuAssertTrue(tc, paf->target_end   == 15);
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 10);
    paf_destruct(paf);
}

static void test_paf_trim_ends_same_strand(CuTest *tc) {
    /* 10M, trim 2 from each end → 6M remains */
    Paf *paf = make_paf("q", 100, 0, 10, true, "t", 100, 0, 10, 10, 10, 60, "10M");
    paf_trim_ends(paf, 2);
    CuAssertTrue(tc, paf->query_start  == 2);
    CuAssertTrue(tc, paf->query_end    == 8);
    CuAssertTrue(tc, paf->target_start == 2);
    CuAssertTrue(tc, paf->target_end   == 8);
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 6);
    paf_destruct(paf);
}

static void test_paf_trim_ends_with_gaps(CuTest *tc) {
    /* "2M1I5M": query span=2+1+5=8, target span=2+5=7
     * Front trim 3: consume 2M (2 aligned) + 1I (gap) + 1 base of 5M → 4M left
     *   query_start += 2+1+1=4, target_start += 2+1=3
     * Back trim 3 on 4M: remove 3 from back → 1M
     *   query_end = 8-3=5, target_end = 7-3=4 */
    Paf *paf = make_paf("q", 100, 0, 8, true, "t", 100, 0, 7, 7, 8, 60, "2M1I5M");
    paf_trim_ends(paf, 3);
    CuAssertTrue(tc, paf->query_start  == 4);
    CuAssertTrue(tc, paf->target_start == 3);
    CuAssertTrue(tc, paf->query_end    == 5);
    CuAssertTrue(tc, paf->target_end   == 4);
    paf_destruct(paf);
}

static void test_paf_trim_end_fraction(CuTest *tc) {
    /* 10M, fraction=0.4 → end_trim = floor(10*0.4/2) = 2 */
    Paf *paf = make_paf("q", 100, 0, 10, true, "t", 100, 0, 10, 10, 10, 60, "10M");
    paf_trim_end_fraction(paf, 0.4f);
    CuAssertTrue(tc, paf->query_start  == 2);
    CuAssertTrue(tc, paf->query_end    == 8);
    CuAssertTrue(tc, paf->target_start == 2);
    CuAssertTrue(tc, paf->target_end   == 8);
    paf_destruct(paf);
}

/* ---- 10. Shatter ---- */

static void test_paf_shatter_single_match(CuTest *tc) {
    Paf *paf = make_paf("q", 100, 0, 5, true, "t", 100, 0, 5, 5, 5, 60, "5M");
    stList *shards = paf_shatter(paf);
    CuAssertIntEquals(tc, 1, stList_length(shards));
    Paf *s = stList_get(shards, 0);
    CuAssertStrEquals(tc, "q", s->query_name);
    CuAssertTrue(tc, s->query_start  == 0);
    CuAssertTrue(tc, s->query_end    == 5);
    CuAssertTrue(tc, s->target_start == 0);
    CuAssertTrue(tc, s->target_end   == 5);
    stList_destruct(shards);
    paf_destruct(paf);
}

static void test_paf_shatter_multi_match(CuTest *tc) {
    /* "3M2D4M": query span=3+4=7, target span=3+2+4=9 */
    Paf *paf = make_paf("q", 100, 0, 7, true, "t", 100, 0, 9, 7, 9, 60, "3M2D4M");
    stList *shards = paf_shatter(paf);
    CuAssertIntEquals(tc, 2, stList_length(shards));

    Paf *s0 = stList_get(shards, 0);
    CuAssertTrue(tc, s0->query_start  == 0);
    CuAssertTrue(tc, s0->query_end    == 3);
    CuAssertTrue(tc, s0->target_start == 0);
    CuAssertTrue(tc, s0->target_end   == 3);

    /* target skips 2D gap: 3+2=5 */
    Paf *s1 = stList_get(shards, 1);
    CuAssertTrue(tc, s1->query_start  == 3);
    CuAssertTrue(tc, s1->query_end    == 7);
    CuAssertTrue(tc, s1->target_start == 5);
    CuAssertTrue(tc, s1->target_end   == 9);

    stList_destruct(shards);
    paf_destruct(paf);
}

static void test_paf_shatter_opposite_strand(CuTest *tc) {
    /* "3M2D4M": query span=7, target span=9; opposite strand.
     * query_coordinate starts at query_end=7 and decrements.
     * 3M: coord 7-3=4 → shard(4,0,3): qs=4,qe=7, ts=0,te=3; target_coord → 3
     * 2D: target_coord → 5
     * 4M: coord 4-4=0 → shard(0,5,4): qs=0,qe=4, ts=5,te=9 */
    Paf *paf = make_paf("q", 100, 0, 7, false, "t", 100, 0, 9, 7, 9, 60, "3M2D4M");
    stList *shards = paf_shatter(paf);
    CuAssertIntEquals(tc, 2, stList_length(shards));

    Paf *s0 = stList_get(shards, 0);
    CuAssertTrue(tc, s0->query_start  == 4);
    CuAssertTrue(tc, s0->query_end    == 7);
    CuAssertTrue(tc, s0->target_start == 0);
    CuAssertTrue(tc, s0->target_end   == 3);

    Paf *s1 = stList_get(shards, 1);
    CuAssertTrue(tc, s1->query_start  == 0);
    CuAssertTrue(tc, s1->query_end    == 4);
    CuAssertTrue(tc, s1->target_start == 5);
    CuAssertTrue(tc, s1->target_end   == 9);

    stList_destruct(shards);
    paf_destruct(paf);
}

/* ---- 11. Mismatch encoding ---- */

static void test_paf_encode_mismatches_all_match(CuTest *tc) {
    Paf *paf = make_paf("q", 5, 0, 5, true, "t", 5, 0, 5, 5, 5, 60, "5M");
    char query_seq[]  = "AAAAA";
    char target_seq[] = "AAAAA";
    paf_encode_mismatches(paf, query_seq, target_seq);
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_match, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    paf_destruct(paf);
}

static void test_paf_encode_mismatches_all_mismatch(CuTest *tc) {
    Paf *paf = make_paf("q", 5, 0, 5, true, "t", 5, 0, 5, 0, 5, 60, "5M");
    char query_seq[]  = "AAAAA";
    char target_seq[] = "CCCCC";
    paf_encode_mismatches(paf, query_seq, target_seq);
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_mismatch, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    paf_destruct(paf);
}

static void test_paf_encode_mismatches_mixed(CuTest *tc) {
    /* target="AACC", query="AATT": AA match, CC vs TT mismatch → 2=2X */
    Paf *paf = make_paf("q", 4, 0, 4, true, "t", 4, 0, 4, 2, 4, 60, "4M");
    char query_seq[]  = "AATT";
    char target_seq[] = "AACC";
    paf_encode_mismatches(paf, query_seq, target_seq);
    CuAssertIntEquals(tc, 2, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_match,    cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 2);
    CuAssertIntEquals(tc, sequence_mismatch, cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 2);
    paf_destruct(paf);
}

static void test_paf_remove_mismatches(CuTest *tc) {
    /* "3=2X1I": = and X merge into 5M; I preserved → "5M1I"
     * query span=3+2+1=6, target span=3+2=5 */
    Paf *paf = make_paf("q", 100, 0, 6, true, "t", 100, 0, 5, 5, 6, 60, "3=2X1I");
    paf_remove_mismatches(paf);
    CuAssertIntEquals(tc, 2, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, match,        cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    CuAssertIntEquals(tc, query_insert, cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 1);
    paf_destruct(paf);
}

/* ---- 12. Coverage tracking ---- */

static void test_coverage_tracking(CuTest *tc) {
    stHash *h = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                  NULL, (void(*)(void*))sequenceCountArray_destruct);

    /* "3M" covering query positions 2,3,4 */
    Paf *paf = make_paf("seq1", 10, 2, 5, true, "t", 100, 0, 3, 3, 3, 60, "3M");

    /* First call creates the array */
    SequenceCountArray *arr1 = get_alignment_count_array(h, paf);
    CuAssertTrue(tc, arr1 != NULL);
    CuAssertTrue(tc, arr1->length == 10);

    /* Second call with same query name returns the same array */
    SequenceCountArray *arr2 = get_alignment_count_array(h, paf);
    CuAssertTrue(tc, arr1 == arr2);

    /* Increment counts */
    increase_alignment_level_counts(arr1, paf);
    CuAssertTrue(tc, arr1->counts[0] == 0);
    CuAssertTrue(tc, arr1->counts[1] == 0);
    CuAssertTrue(tc, arr1->counts[2] == 1);
    CuAssertTrue(tc, arr1->counts[3] == 1);
    CuAssertTrue(tc, arr1->counts[4] == 1);
    CuAssertTrue(tc, arr1->counts[5] == 0);

    paf_destruct(paf);
    stHash_destruct(h);
}

/* ---- 13. Interval functions ---- */

static void test_decode_fasta_header(CuTest *tc) {
    /* fasta_chunk format: "name|sequenceLength|chunkStart"
     * decode pops last two fields as start then length */
    Interval *iv = decode_fasta_header("seqname|100|0");
    CuAssertStrEquals(tc, "seqname", iv->name);
    CuAssertTrue(tc, iv->start  == 0);
    CuAssertTrue(tc, iv->length == 100);
    interval_destruct(iv);
}

static void test_cmp_intervals(CuTest *tc) {
    Interval a = { .name = "chr1", .start = 10, .end = 100, .length = 90 };
    Interval b = { .name = "chr1", .start = 20, .end = 200, .length = 180 };
    Interval c = { .name = "chr2", .start = 5,  .end = 50,  .length = 45 };

    /* same name: compare by start */
    CuAssertTrue(tc, cmp_intervals(&a, &b) < 0);
    CuAssertTrue(tc, cmp_intervals(&b, &a) > 0);
    CuAssertTrue(tc, cmp_intervals(&a, &a) == 0);

    /* different names: lexicographic */
    CuAssertTrue(tc, cmp_intervals(&a, &c) < 0); /* "chr1" < "chr2" */
    CuAssertTrue(tc, cmp_intervals(&c, &a) > 0);
}

/* ---- 14. paf_trim_unreliable_tails ---- */

static void test_paf_trim_unreliable_tails_trims_tails(CuTest *tc) {
    /* "2X5=2X": query span=9, target span=9, same_strand=false.
     * matches=5, mismatches=4, identity=5/9; score_fraction=0 so threshold=identity.
     * For same_strand=false the suffix is trimmed via invert+trim+invert:
     * both 2X tails are removed, leaving "5=", query=(2,7), target=(2,7). */
    Paf *paf = make_paf("q", 9, 0, 9, false, "t", 9, 0, 9, 5, 9, 60, "2X5=2X");
    paf_trim_unreliable_tails(paf, 0.0f, 1.0f);
    CuAssertTrue(tc, paf->query_start  == 2);
    CuAssertTrue(tc, paf->query_end    == 7);
    CuAssertTrue(tc, paf->target_start == 2);
    CuAssertTrue(tc, paf->target_end   == 7);
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_match, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    paf_destruct(paf);
}

static void test_paf_trim_unreliable_tails_no_trim(CuTest *tc) {
    /* score_fraction=1.0 → identity_threshold = identity - identity*1 = 0.
     * No prefix can have identity < 0, so nothing is trimmed. */
    Paf *paf = make_paf("q", 9, 0, 9, true, "t", 9, 0, 9, 5, 9, 60, "2X5=2X");
    paf_trim_unreliable_tails(paf, 1.0f, 1.0f);
    CuAssertTrue(tc, paf->query_start  == 0);
    CuAssertTrue(tc, paf->query_end    == 9);
    CuAssertTrue(tc, paf->target_start == 0);
    CuAssertTrue(tc, paf->target_end   == 9);
    CuAssertIntEquals(tc, 3, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_mismatch, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 2);
    CuAssertIntEquals(tc, sequence_match, cigar_get(paf->cigar, 1)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 1)->length == 5);
    CuAssertIntEquals(tc, sequence_mismatch, cigar_get(paf->cigar, 2)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 2)->length == 2);
    paf_destruct(paf);
}

static void test_paf_trim_unreliable_tails_opposite_strand(CuTest *tc) {
    /* same_strand=false: paf_trim_upto decrements query_end (not query_start)
     * because the query is on the opposite strand.
     * "2X5=": query span=7, target span=7, same_strand=false.
     * identity=5/7; threshold=5/7 (score_fraction=0).
     * After prefix trim: query_end decremented by 2 to 5; target_start incremented to 2.
     * The suffix "5=" has identity 1.0 > threshold and is not trimmed. */
    Paf *paf = make_paf("q", 9, 0, 7, false, "t", 9, 0, 7, 5, 7, 60, "2X5=");
    paf_trim_unreliable_tails(paf, 0.0f, 1.0f);
    CuAssertTrue(tc, paf->query_start  == 0);  /* not incremented for opp strand */
    CuAssertTrue(tc, paf->query_end    == 5);  /* decremented by 2 */
    CuAssertTrue(tc, paf->target_start == 2);  /* incremented */
    CuAssertTrue(tc, paf->target_end   == 7);  /* unchanged */
    CuAssertIntEquals(tc, 1, cigar_count(paf->cigar));
    CuAssertIntEquals(tc, sequence_match, cigar_get(paf->cigar, 0)->op);
    CuAssertTrue(tc, cigar_get(paf->cigar, 0)->length == 5);
    paf_destruct(paf);
}

/* ---- 15. paf_pretty_print ---- */

static void test_paf_pretty_print_basic(CuTest *tc) {
    /* Call paf_pretty_print with include_alignment=false and verify that
     * a non-empty header line is written to the output file. */
    Paf *paf = make_paf("q", 10, 0, 5, true, "t", 10, 0, 5, 5, 5, 60, "5=");
    FILE *fh = tmpfile();
    CuAssertTrue(tc, fh != NULL);
    paf_pretty_print(paf, NULL, NULL, fh, false);
    CuAssertTrue(tc, ftell(fh) > 0);  /* output is non-empty */
    fclose(fh);
    paf_destruct(paf);
}

/* ---- 16. paf_check (positive path) ---- */

static void test_paf_check_valid(CuTest *tc) {
    /* paf_check should not abort on valid records; reaching the end of this
     * function without st_errAbort is the success criterion. */
    Paf *p;

    /* same-strand, no cigar */
    p = make_paf("q", 100, 0, 50, true, "t", 200, 10, 60, 50, 50, 60, NULL);
    paf_check(p);
    paf_destruct(p);

    /* opposite-strand, no cigar */
    p = make_paf("q", 100, 0, 50, false, "t", 200, 10, 60, 50, 50, 60, NULL);
    paf_check(p);
    paf_destruct(p);

    /* same-strand with simple cigar */
    p = make_paf("q", 100, 0, 5, true, "t", 100, 0, 5, 5, 5, 60, "5=");
    paf_check(p);
    paf_destruct(p);

    /* cigar with indels: "3=2X1I2D", query span=3+2+1=6, target span=3+2+2=7 */
    p = make_paf("q", 100, 0, 6, true, "t", 100, 0, 7, 5, 8, 60, "3=2X1I2D");
    paf_check(p);
    paf_destruct(p);

    CuAssertTrue(tc, 1);  /* reached here without aborting */
}

/* ---- 17. Left-aligning gaps ---- */

static char *cigar_to_str(Cigar *c, bool reversed) {
    static const char ops[] = { 'M', 'I', 'D', '=', 'X' };
    stList *parts = stList_construct3(0, free);
    for (int64_t i = 0; i < cigar_count(c); i++) {
        CigarRecord *r = cigar_get(c, reversed ? cigar_count(c) - 1 - i : i);
        stList_append(parts, stString_print("%" PRIi64 "%c", (int64_t)r->length, ops[r->op]));
    }
    char *s = stString_join2("", parts);
    stList_destruct(parts);
    return s;
}

static char *reverse_str(const char *s) {
    int64_t n = strlen(s);
    char *r = st_malloc(n + 1);
    for (int64_t i = 0; i < n; i++) {
        r[i] = s[n - 1 - i];
    }
    r[n] = '\0';
    return r;
}

/* The target and the query as aligned get flanks, so that the alignment does not start at zero, and the query is
 * reverse complemented on the opposite strand */
static Paf *make_flanked_paf(const char *target, const char *aligned_query, bool same_strand, const char *cigar,
                             char **t, char **q) {
    *t = stString_print("TTTT%sTTT", target);
    char *full_q = stString_print("CCC%sGG", aligned_query);
    *q = same_strand ? stString_copy(full_q) : stString_reverseComplementString(full_q);
    free(full_q);
    int64_t ql = strlen(aligned_query), tl = strlen(target), qs = same_strand ? 3 : 2;
    return make_paf("q", ql + 5, qs, qs + ql, same_strand, "t", tl + 7, 4, 4 + tl, 0, 0, 60, cigar);
}

/* Counts of the alignment's columns of each kind (as paf_stats_calc), which are what its score is made of */
static void alignment_stats(const char *target, const char *aligned_query, bool same_strand, const char *cigar,
                            int64_t *stats) {
    char *t, *q;
    Paf *p = make_flanked_paf(target, aligned_query, same_strand, cigar, &t, &q);
    paf_check(p);
    paf_remove_mismatches(p);
    paf_encode_mismatches(p, q, t);
    paf_stats_calc(p, &stats[0], &stats[1], &stats[2], &stats[3], &stats[4], &stats[5], 1);
    paf_destruct(p);
    free(t); free(q);
}

/* Left-align the cigar, returning the new one, and the number of gaps moved in moved if not NULL */
static char *gap_align_str(const char *target, const char *aligned_query, bool same_strand, const char *cigar,
                           int64_t *moved, bool canonical) {
    char *t, *q;
    Paf *p = make_flanked_paf(target, aligned_query, same_strand, cigar, &t, &q);
    int64_t m = canonical ? paf_canonical_align(p, q, t) : paf_left_align(p, q, t);
    if (moved) *moved = m;
    paf_check(p);
    char *s = cigar_to_str(p->cigar, 0);
    paf_destruct(p);
    free(t); free(q);
    return s;
}

static char *left_align_str(const char *target, const char *aligned_query, bool same_strand, const char *cigar,
                            int64_t *moved) {
    return gap_align_str(target, aligned_query, same_strand, cigar, moved, 0);
}

static char *canonical_align_str(const char *target, const char *aligned_query, bool same_strand, const char *cigar) {
    return gap_align_str(target, aligned_query, same_strand, cigar, NULL, 1);
}

static void check_left_align(CuTest *tc, const char *target, const char *aligned_query, const char *cigar,
                             const char *expected) {
    for (int strand = 0; strand < 2; strand++) {
        char *got = left_align_str(target, aligned_query, strand == 0, cigar, NULL);
        CuAssertStrEquals(tc, expected, got);
        free(got);
    }
}

static void test_paf_left_align_tandem_delete(CuTest *tc) {
    /* one CA of (CA)3 deleted, written at the right end of the repeat */
    check_left_align(tc, "TTGCACACAGTT", "TTGCACAGTT", "7M2D3M", "3M2D7M");
    /* already left: nothing moves */
    check_left_align(tc, "TTGCACACAGTT", "TTGCACAGTT", "3M2D7M", "3M2D7M");
}

static void test_paf_left_align_tandem_insert(CuTest *tc) {
    check_left_align(tc, "TTGCACAGTT", "TTGCACACAGTT", "7M2I3M", "3M2I7M");
}

static void test_paf_left_align_homopolymer_encoded(CuTest *tc) {
    check_left_align(tc, "GCAAAATC", "GCAAATC", "5=1D2=", "2=1D5=");
    /* a mismatch before the repeat stops nothing it shouldn't */
    check_left_align(tc, "GCAAAATC", "GTAAATC", "1=1X3=1D2=", "1=1X1D5=");
    /* a mismatch inside the repeat moves to the other side of the gap, keeping its pair of bases */
    check_left_align(tc, "GAAAAC", "GATAC", "2=1X1=1D1=", "1=1D1=1X2=");
    /* M in, M out */
    check_left_align(tc, "GAAAAC", "GATAC", "4M1D1M", "1M1D4M");
}

static void test_paf_left_align_alignment_start(CuTest *tc) {
    /* the gap could go to the very start, but stops one column short of it */
    check_left_align(tc, "AAAT", "AAT", "2M1D1M", "1M1D2M");
    check_left_align(tc, "AAT", "AAAT", "2M1I1M", "1M1I2M");
}

static void test_paf_left_align_merge(CuTest *tc) {
    /* two single-base deletes in a run of A meet, and become one */
    check_left_align(tc, "GAAAAT", "GAAT", "2M1D1M1D1M", "1M2D3M");
    /* adjacent gaps of the same kind are merged even when nothing moves */
    check_left_align(tc, "GCTA", "GA", "1M1D1D1M", "1M2D1M");
    /* a gap does not run into one of the other kind */
    check_left_align(tc, "GAAT", "GCAT", "1M1I1D2M", "1M1I1D2M");
}

/* a small deterministic generator, so the random tests do not depend on, or disturb, the global one */
static uint64_t la_state = 88172645463325252ULL;
static int64_t la_rand(int64_t n) {
    la_state ^= la_state << 13; la_state ^= la_state >> 7; la_state ^= la_state << 17;
    return (int64_t)(la_state % (uint64_t)n);
}

static char la_base(void) { // mostly A and C, so there are plenty of repeats for gaps to slide in
    return la_rand(8) == 0 ? "GT"[la_rand(2)] : "AC"[la_rand(2)];
}

static void test_paf_left_align_equivalent_placements(CuTest *tc) {
    /* Any of the places a gap can equivalently go must left-align to the same one. A right-aligned placement is
     * made by left-aligning the reversed alignment, and must come back to the left-aligned one. */
    for (int64_t test = 0; test < 2000; test++) {
        int64_t n = 10 + la_rand(30), l = 1 + la_rand(4), p = 1 + la_rand(n - 1 - l);
        bool del = la_rand(2);
        char *t = st_malloc(n + 1);
        for (int64_t i = 0; i < n; i++) t[i] = la_base();
        t[n] = '\0';
        char *q = st_malloc(n + l + 1);
        int64_t ql = 0;
        for (int64_t i = 0; i < n; i++) {
            if (!del && i == p) for (int64_t j = 0; j < l; j++) q[ql++] = la_base();
            if (!del || i < p || i >= p + l) q[ql++] = t[i];
        }
        q[ql] = '\0';
        for (int64_t k = la_rand(3); k > 0; k--) { // a few substitutions, away from inserted bases
            int64_t i = la_rand(ql);
            if (del || i < p || i >= p + l) q[i] = la_base();
        }
        char *cigar = del ? stString_print("%" PRIi64 "M%" PRIi64 "D%" PRIi64 "M", p, l, n - p - l)
                          : stString_print("%" PRIi64 "M%" PRIi64 "I%" PRIi64 "M", p, l, n - p);
        bool same_strand = la_rand(2);

        int64_t stats0[6], stats1[6], moved;
        char *left = left_align_str(t, q, same_strand, cigar, NULL);
        char *rt = reverse_str(t), *rq = reverse_str(q);
        Cigar *rc = cigar_parse(cigar);
        char *rcigar = cigar_to_str(rc, 1);
        char *rleft = left_align_str(rt, rq, same_strand, rcigar, NULL);
        Cigar *rlc = cigar_parse(rleft);
        char *right = cigar_to_str(rlc, 1);
        char *left2 = left_align_str(t, q, same_strand, right, NULL);
        CuAssertStrEquals(tc, left, left2);
        alignment_stats(t, q, same_strand, cigar, stats0);
        alignment_stats(t, q, same_strand, right, stats1);
        for (int i = 0; i < 6; i++) CuAssertIntEquals(tc, (int)stats0[i], (int)stats1[i]);
        char *left3 = left_align_str(t, q, same_strand, left, &moved);
        CuAssertStrEquals(tc, left, left3);
        CuAssertIntEquals(tc, 0, (int)moved);

        free(t); free(q); free(cigar); free(left); free(rt); free(rq); cigar_destruct(rc); free(rcigar);
        free(rleft); cigar_destruct(rlc); free(right); free(left2); free(left3);
    }
}

static void test_paf_left_align_preserves_columns(CuTest *tc) {
    /* Alignments with many gaps: whatever moves, the counts of matches, mismatches and gap bases are unchanged,
     * and a second pass moves nothing */
    for (int64_t test = 0; test < 2000; test++) {
        stList *ops = stList_construct3(0, free);
        char t[512], q[512];
        int64_t tl = 0, ql = 0;
        int64_t m = 1 + la_rand(4);
        for (int64_t i = 0; i < m; i++) { t[tl] = q[ql] = la_base(); tl++; ql++; }
        stList_append(ops, stString_print("%" PRIi64 "M", m));
        for (int64_t k = 1 + la_rand(6); k > 0; k--) {
            int64_t l = 1 + la_rand(3);
            if (la_rand(2)) { for (int64_t i = 0; i < l; i++) t[tl++] = la_base(); stList_append(ops, stString_print("%" PRIi64 "D", l)); }
            else { for (int64_t i = 0; i < l; i++) q[ql++] = la_base(); stList_append(ops, stString_print("%" PRIi64 "I", l)); }
            m = la_rand(5); // no run between gaps sometimes
            for (int64_t i = 0; i < m; i++) { t[tl] = la_base(); q[ql] = la_rand(5) ? t[tl] : la_base(); tl++; ql++; }
            if (m > 0) stList_append(ops, stString_print("%" PRIi64 "M", m));
        }
        t[tl] = q[ql] = '\0';
        char *cigar = stString_join2("", ops);
        bool same_strand = la_rand(2);

        int64_t before[6], after[6], moved;
        alignment_stats(t, q, same_strand, cigar, before);
        char *left = left_align_str(t, q, same_strand, cigar, NULL);
        alignment_stats(t, q, same_strand, left, after);
        CuAssertIntEquals(tc, (int)before[0], (int)after[0]); // matches
        CuAssertIntEquals(tc, (int)before[1], (int)after[1]); // mismatches
        CuAssertIntEquals(tc, (int)before[4], (int)after[4]); // inserted bases
        CuAssertIntEquals(tc, (int)before[5], (int)after[5]); // deleted bases
        CuAssertTrue(tc, after[2] + after[3] <= before[2] + before[3]); // gaps only ever merge
        char *left2 = left_align_str(t, q, same_strand, left, &moved);
        CuAssertStrEquals(tc, left, left2);
        CuAssertIntEquals(tc, 0, (int)moved);

        // the same of canonical placement, which a second pass leaves where it is (it moves gaps left and back)
        char *canon = canonical_align_str(t, q, same_strand, cigar);
        alignment_stats(t, q, same_strand, canon, after);
        CuAssertIntEquals(tc, (int)before[0], (int)after[0]);
        CuAssertIntEquals(tc, (int)before[1], (int)after[1]);
        CuAssertIntEquals(tc, (int)before[4], (int)after[4]);
        CuAssertIntEquals(tc, (int)before[5], (int)after[5]);
        char *canon2 = canonical_align_str(t, q, same_strand, canon);
        CuAssertStrEquals(tc, canon, canon2);

        free(cigar); free(left); free(left2); free(canon); free(canon2);
        stList_destruct(ops);
    }
}


static void test_paf_canonical_align_tandem(CuTest *tc) {
    for (int strand = 0; strand < 2; strand++) {
        char *c;
        /* (CA)3 reads smaller than (TG)3, so a deletion in it stays at the left of CACACA ... */
        c = canonical_align_str("TTGCACACAGTT", "TTGCACAGTT", strand, "7M2D3M");
        CuAssertStrEquals(tc, "3M2D7M", c); free(c);
        /* ... and seen from the other strand, as (TG)3, it goes right: the same bases either way */
        c = canonical_align_str("AACTGTGTGCAA", "AACTGTGCAA", strand, "3M2D7M");
        CuAssertStrEquals(tc, "7M2D3M", c); free(c);
        /* an insert is judged on the query */
        c = canonical_align_str("AACTGTGCAA", "AACTGTGTGCAA", strand, "3M2I7M");
        CuAssertStrEquals(tc, "7M2I3M", c); free(c);
        /* (AT)3 is its own reverse complement, so the flanks decide: TTG..GTT against its reverse complement
         * AAC..CAA is bigger, and the gap goes right */
        c = canonical_align_str("TTGATATATGTT", "TTGATATGTT", strand, "3M2D7M");
        CuAssertStrEquals(tc, "7M2D3M", c); free(c);
    }
}

static void test_paf_canonical_align_strand_invariant(CuTest *tc) {
    /* An indel must land on the same bases seen from either strand: aligning the reverse complements of both
     * sequences, with the cigar reversed, must give the reversed cigar. With left_align it would not.
     * The one exception is a gap whose surroundings are their own reverse complement out to the ends of the
     * alignment, which looks the same from both strands and so cannot be told apart: both stay left */
    int64_t left_differs = 0, unsettled = 0;
    for (int64_t test = 0; test < 2000; test++) {
        int64_t n = 10 + la_rand(30), l = 1 + la_rand(4), p = 1 + la_rand(n - 1 - l);
        bool del = la_rand(2);
        char *t = st_malloc(n + 1);
        for (int64_t i = 0; i < n; i++) t[i] = "ACGT"[la_rand(4)];
        t[n] = '\0';
        char *q = st_malloc(n + l + 1);
        int64_t ql = 0;
        for (int64_t i = 0; i < n; i++) {
            if (!del && i == p) for (int64_t j = 0; j < l; j++) q[ql++] = t[p - l + j >= 0 ? p - l + j : 0]; // a duplication
            if (!del || i < p || i >= p + l) q[ql++] = t[i];
        }
        q[ql] = '\0';
        char *cigar = del ? stString_print("%" PRIi64 "M%" PRIi64 "D%" PRIi64 "M", p, l, n - p - l)
                          : stString_print("%" PRIi64 "M%" PRIi64 "I%" PRIi64 "M", p, l, n - p);
        bool same_strand = la_rand(2);
        char *rt = stString_reverseComplementString(t), *rq = stString_reverseComplementString(q);
        Cigar *c = cigar_parse(cigar);
        char *rcigar = cigar_to_str(c, 1);

        char *canon = canonical_align_str(t, q, same_strand, cigar);
        char *rcanon = canonical_align_str(rt, rq, same_strand, rcigar);
        Cigar *rc = cigar_parse(rcanon);
        char *rcanon_back = cigar_to_str(rc, 1);

        char *left = left_align_str(t, q, same_strand, cigar, NULL);
        char *rleft = left_align_str(rt, rq, same_strand, rcigar, NULL);
        Cigar *rl = cigar_parse(rleft);
        char *rleft_back = cigar_to_str(rl, 1);
        left_differs += strcmp(left, rleft_back) != 0;

        if (strcmp(canon, left) == 0 && strcmp(rcanon, rleft) == 0 && strcmp(left, rleft_back) != 0) {
            unsettled++; // left from both strands, which is only right when the strands cannot be told apart
        } else {
            CuAssertStrEquals(tc, canon, rcanon_back);
        }

        free(t); free(q); free(cigar); free(rt); free(rq); cigar_destruct(c); free(rcigar); free(canon);
        free(rcanon); cigar_destruct(rc); free(rcanon_back); free(left); free(rleft); cigar_destruct(rl);
        free(rleft_back);
    }
    CuAssertTrue(tc, left_differs > 100); // the test has teeth: plain left-aligning is not strand invariant
    CuAssertTrue(tc, unsettled < left_differs / 20);
}

/* ---- Registration ---- */

CuSuite *addPafUnitTestSuite(void) {
    CuSuite *suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_cigar_parse_empty);
    SUITE_ADD_TEST(suite, test_cigar_parse_single);
    SUITE_ADD_TEST(suite, test_cigar_parse_all_ops);
    SUITE_ADD_TEST(suite, test_cigar_parse_large_length);
    SUITE_ADD_TEST(suite, test_cigar_count_get);
    SUITE_ADD_TEST(suite, test_paf_parse_minimal);
    SUITE_ADD_TEST(suite, test_paf_parse_with_cigar);
    SUITE_ADD_TEST(suite, test_paf_parse_cigar_string_mode);
    SUITE_ADD_TEST(suite, test_paf_parse_optional_tags);
    SUITE_ADD_TEST(suite, test_paf_parse_strand);
    SUITE_ADD_TEST(suite, test_paf_roundtrip_no_cigar);
    SUITE_ADD_TEST(suite, test_paf_roundtrip_with_cigar);
    SUITE_ADD_TEST(suite, test_paf_read_write);
    SUITE_ADD_TEST(suite, test_read_write_pafs_list);
    SUITE_ADD_TEST(suite, test_paf_stats_calc_all_match);
    SUITE_ADD_TEST(suite, test_paf_stats_calc_mixed);
    SUITE_ADD_TEST(suite, test_paf_stats_calc_zero_flag);
    SUITE_ADD_TEST(suite, test_paf_invert_same_strand);
    SUITE_ADD_TEST(suite, test_paf_invert_opposite_strand);
    SUITE_ADD_TEST(suite, test_paf_invert_double);
    SUITE_ADD_TEST(suite, test_aligned_bases);
    SUITE_ADD_TEST(suite, test_paf_trim_ends_zero);
    SUITE_ADD_TEST(suite, test_paf_trim_ends_same_strand);
    SUITE_ADD_TEST(suite, test_paf_trim_ends_with_gaps);
    SUITE_ADD_TEST(suite, test_paf_trim_end_fraction);
    SUITE_ADD_TEST(suite, test_paf_shatter_single_match);
    SUITE_ADD_TEST(suite, test_paf_shatter_multi_match);
    SUITE_ADD_TEST(suite, test_paf_shatter_opposite_strand);
    SUITE_ADD_TEST(suite, test_paf_encode_mismatches_all_match);
    SUITE_ADD_TEST(suite, test_paf_encode_mismatches_all_mismatch);
    SUITE_ADD_TEST(suite, test_paf_encode_mismatches_mixed);
    SUITE_ADD_TEST(suite, test_paf_remove_mismatches);
    SUITE_ADD_TEST(suite, test_coverage_tracking);
    SUITE_ADD_TEST(suite, test_decode_fasta_header);
    SUITE_ADD_TEST(suite, test_cmp_intervals);
    SUITE_ADD_TEST(suite, test_paf_trim_unreliable_tails_trims_tails);
    SUITE_ADD_TEST(suite, test_paf_trim_unreliable_tails_no_trim);
    SUITE_ADD_TEST(suite, test_paf_trim_unreliable_tails_opposite_strand);
    SUITE_ADD_TEST(suite, test_paf_pretty_print_basic);
    SUITE_ADD_TEST(suite, test_paf_check_valid);
    SUITE_ADD_TEST(suite, test_paf_left_align_tandem_delete);
    SUITE_ADD_TEST(suite, test_paf_left_align_tandem_insert);
    SUITE_ADD_TEST(suite, test_paf_left_align_homopolymer_encoded);
    SUITE_ADD_TEST(suite, test_paf_left_align_alignment_start);
    SUITE_ADD_TEST(suite, test_paf_left_align_merge);
    SUITE_ADD_TEST(suite, test_paf_left_align_equivalent_placements);
    SUITE_ADD_TEST(suite, test_paf_left_align_preserves_columns);
    SUITE_ADD_TEST(suite, test_paf_canonical_align_tandem);
    SUITE_ADD_TEST(suite, test_paf_canonical_align_strand_invariant);
    return suite;
}
