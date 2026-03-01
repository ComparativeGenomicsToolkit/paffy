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
    return suite;
}
