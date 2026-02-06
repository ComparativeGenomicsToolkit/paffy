#include "paf.h"
#include "CuTest.h"
#include "sonLib.h"

/*
 * Unit tests for paffy split_file
 */

static void write_test_paf_file(const char *path) {
    FILE *fh = fopen(path, "w");
    assert(fh != NULL);
    fprintf(fh, "q1\t100\t0\t50\t+\tchr1\t1000\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q2\t100\t0\t50\t+\tchr2\t500\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q3\t100\t0\t50\t-\tchr1\t1000\t100\t150\t50\t50\t60\n");
    fprintf(fh, "q4\t100\t0\t50\t+\tsmall_a\t300\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q5\t100\t0\t50\t+\tsmall_b\t200\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q6\t100\t0\t50\t+\tsmall_c\t400\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q7\t100\t0\t50\t+\tsmall_a\t300\t50\t100\t50\t50\t60\n");
    fprintf(fh, "q8\t100\t0\t50\t+\tsmall_d\t150\t0\t50\t50\t50\t60\n");
    fclose(fh);
}

static int file_exists(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh != NULL) {
        fclose(fh);
        return 1;
    }
    return 0;
}

static int count_records(const char *path) {
    FILE *fh = fopen(path, "r");
    if (fh == NULL) return -1;
    stList *pafs = read_pafs(fh, 0);
    int count = stList_length(pafs);
    fclose(fh);
    stList_destruct(pafs);
    return count;
}

static void test_split_file_basic(CuTest *testCase) {
    // Write test PAF
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_";
    write_test_paf_file(input);

    // Run split_file with no minTargetLength
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s", input, prefix) == 0);

    // Verify all expected files exist
    CuAssertTrue(testCase, file_exists("./tests/temp_split_chr1.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_chr2.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_a.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_b.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_c.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_d.paf"));

    // Verify record counts
    CuAssertIntEquals(testCase, 2, count_records("./tests/temp_split_chr1.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_chr2.paf"));
    CuAssertIntEquals(testCase, 2, count_records("./tests/temp_split_small_a.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_b.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_c.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_d.paf"));

    // Verify target names in chr1 file
    FILE *fh = fopen("./tests/temp_split_chr1.paf", "r");
    stList *pafs = read_pafs(fh, 0);
    fclose(fh);
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        CuAssertTrue(testCase, strcmp(paf->target_name, "chr1") == 0);
    }
    stList_destruct(pafs);

    // Verify target names in small_a file
    fh = fopen("./tests/temp_split_small_a.paf", "r");
    pafs = read_pafs(fh, 0);
    fclose(fh);
    CuAssertIntEquals(testCase, 2, stList_length(pafs));
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        CuAssertTrue(testCase, strcmp(paf->target_name, "small_a") == 0);
    }
    stList_destruct(pafs);

    // Total records across all files
    int total = count_records("./tests/temp_split_chr1.paf") +
                count_records("./tests/temp_split_chr2.paf") +
                count_records("./tests/temp_split_small_a.paf") +
                count_records("./tests/temp_split_small_b.paf") +
                count_records("./tests/temp_split_small_c.paf") +
                count_records("./tests/temp_split_small_d.paf");
    CuAssertIntEquals(testCase, 8, total);

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf ./tests/temp_split_chr1.paf "
              "./tests/temp_split_chr2.paf ./tests/temp_split_small_a.paf "
              "./tests/temp_split_small_b.paf ./tests/temp_split_small_c.paf "
              "./tests/temp_split_small_d.paf");
}

static void test_split_file_min_target_length(CuTest *testCase) {
    // Write test PAF
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_";
    write_test_paf_file(input);

    // Run split_file with -m 500
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s -m 500", input, prefix) == 0);

    // Large contigs get their own files
    CuAssertTrue(testCase, file_exists("./tests/temp_split_chr1.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_chr2.paf"));
    CuAssertIntEquals(testCase, 2, count_records("./tests/temp_split_chr1.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_chr2.paf"));

    // Small contigs are grouped into small_N files
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_0.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_1.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_2.paf"));

    // Per-target files should NOT exist for small contigs
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_small_a.paf"));
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_small_b.paf"));
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_small_c.paf"));
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_small_d.paf"));

    // Verify small_0 has small_a (300) + small_b (200) = 3 records (q4, q5, q7)
    CuAssertIntEquals(testCase, 3, count_records("./tests/temp_split_small_0.paf"));

    // Verify small_1 has small_c (400) = 1 record (q6)
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_1.paf"));

    // Verify small_2 has small_d (150) = 1 record (q8)
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_2.paf"));

    // Verify all records in small_0 have target_name in {small_a, small_b}
    FILE *fh = fopen("./tests/temp_split_small_0.paf", "r");
    stList *pafs = read_pafs(fh, 0);
    fclose(fh);
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        CuAssertTrue(testCase, strcmp(paf->target_name, "small_a") == 0 ||
                               strcmp(paf->target_name, "small_b") == 0);
    }
    stList_destruct(pafs);

    // Verify both small_a records (q4, q7) are in small_0 (contig locality preserved)
    fh = fopen("./tests/temp_split_small_0.paf", "r");
    pafs = read_pafs(fh, 0);
    fclose(fh);
    int small_a_count = 0;
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        if (strcmp(paf->target_name, "small_a") == 0) {
            small_a_count++;
        }
    }
    CuAssertIntEquals(testCase, 2, small_a_count);
    stList_destruct(pafs);

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf ./tests/temp_split_chr1.paf "
              "./tests/temp_split_chr2.paf ./tests/temp_split_small_0.paf "
              "./tests/temp_split_small_1.paf ./tests/temp_split_small_2.paf");
}

static void test_split_file_all_small(CuTest *testCase) {
    // Write a PAF with only small contigs
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_";
    FILE *fh = fopen(input, "w");
    assert(fh != NULL);
    fprintf(fh, "q1\t100\t0\t50\t+\tctg1\t100\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q2\t100\t0\t50\t+\tctg2\t100\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q3\t100\t0\t50\t+\tctg3\t100\t0\t50\t50\t50\t60\n");
    fclose(fh);

    // Run with -m 250: ctg1(100)+ctg2(100)=200<=250 -> small_0, ctg3(100): 200+100=300>250 -> small_1
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s -m 250", input, prefix) == 0);

    // Should have small_0 and small_1
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_0.paf"));
    CuAssertTrue(testCase, file_exists("./tests/temp_split_small_1.paf"));

    CuAssertIntEquals(testCase, 2, count_records("./tests/temp_split_small_0.paf"));
    CuAssertIntEquals(testCase, 1, count_records("./tests/temp_split_small_1.paf"));

    // No per-target files should exist
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_ctg1.paf"));
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_ctg2.paf"));
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_ctg3.paf"));

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf ./tests/temp_split_small_0.paf "
              "./tests/temp_split_small_1.paf");
}

static void test_split_file_empty_input(CuTest *testCase) {
    // Write an empty file
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_empty_";
    FILE *fh = fopen(input, "w");
    assert(fh != NULL);
    fclose(fh);

    // Run split_file
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s", input, prefix) == 0);

    // No output files should be created - verify a few likely names don't exist
    CuAssertTrue(testCase, !file_exists("./tests/temp_split_empty_small_0.paf"));

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf");
}

static void test_split_file_single_target(CuTest *testCase) {
    // Write PAF with 3 records all targeting chrX
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_";
    FILE *fh = fopen(input, "w");
    assert(fh != NULL);
    fprintf(fh, "q1\t100\t0\t50\t+\tchrX\t5000\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q2\t100\t0\t50\t+\tchrX\t5000\t100\t150\t50\t50\t60\n");
    fprintf(fh, "q3\t100\t0\t50\t-\tchrX\t5000\t200\t250\t50\t50\t60\n");
    fclose(fh);

    // Run split_file
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s", input, prefix) == 0);

    // Exactly one output file
    CuAssertTrue(testCase, file_exists("./tests/temp_split_chrX.paf"));
    CuAssertIntEquals(testCase, 3, count_records("./tests/temp_split_chrX.paf"));

    // Verify all records have target_name chrX
    fh = fopen("./tests/temp_split_chrX.paf", "r");
    stList *pafs = read_pafs(fh, 0);
    fclose(fh);
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        CuAssertTrue(testCase, strcmp(paf->target_name, "chrX") == 0);
    }
    stList_destruct(pafs);

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf ./tests/temp_split_chrX.paf");
}

static void test_split_file_sanitize_filename(CuTest *testCase) {
    // Write PAF with a target name containing '/'
    const char *input = "./tests/temp_split_input.paf";
    const char *prefix = "./tests/temp_split_";
    FILE *fh = fopen(input, "w");
    assert(fh != NULL);
    fprintf(fh, "q1\t100\t0\t50\t+\tcontig/scaffold_1\t2000\t0\t50\t50\t50\t60\n");
    fprintf(fh, "q2\t100\t0\t50\t+\tcontig/scaffold_1\t2000\t100\t150\t50\t50\t60\n");
    fclose(fh);

    // Run split_file
    CuAssertTrue(testCase, st_system("./bin/paffy split_file -i %s -p %s", input, prefix) == 0);

    // Slashes should be replaced with underscores in filename
    CuAssertTrue(testCase, file_exists("./tests/temp_split_contig_scaffold_1.paf"));
    CuAssertIntEquals(testCase, 2, count_records("./tests/temp_split_contig_scaffold_1.paf"));

    // Verify the original target_name is preserved in the records
    fh = fopen("./tests/temp_split_contig_scaffold_1.paf", "r");
    stList *pafs = read_pafs(fh, 0);
    fclose(fh);
    for (int64_t i = 0; i < stList_length(pafs); i++) {
        Paf *paf = stList_get(pafs, i);
        CuAssertTrue(testCase, strcmp(paf->target_name, "contig/scaffold_1") == 0);
    }
    stList_destruct(pafs);

    // Cleanup
    st_system("rm -f ./tests/temp_split_input.paf ./tests/temp_split_contig_scaffold_1.paf");
}

CuSuite* addPafSplitFileTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_split_file_basic);
    SUITE_ADD_TEST(suite, test_split_file_min_target_length);
    SUITE_ADD_TEST(suite, test_split_file_all_small);
    SUITE_ADD_TEST(suite, test_split_file_empty_input);
    SUITE_ADD_TEST(suite, test_split_file_single_target);
    SUITE_ADD_TEST(suite, test_split_file_sanitize_filename);
    return suite;
}
