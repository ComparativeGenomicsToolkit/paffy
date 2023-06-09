#include "CuTest.h"
#include "sonLib.h"
#include "bioioC.h"

/*
 * Test fasta extract
 */

static char *test_fa_file = "./tests/temp.fa";
static char *test_fasta_chunks_dir = "./tests/temp_fastas";
static char *test_chunks_file = "./tests/chunks.txt";
static char *test_dechunked_fa_file = "./tests/temp2.fa";

static void test_fasta_chunk_and_merge(CuTest *testCase) {
    int64_t chunk_size = 1000000;
    int64_t overlap = 10000;

    // Download a test sequence
    CuAssertTrue(testCase, st_system("wget https://glennhickey.s3.amazonaws.com/share/hg38_preprocessed_chr10.fa -O %s", test_fa_file) == 0);

    // Run fasta chunk to get a list of chunks
    CuAssertTrue(testCase, st_system("faffy chunk --logLevel DEBUG %s -d %s -c %" PRIi64 " -o %" PRIi64 " > %s", test_fa_file, test_fasta_chunks_dir,
                                      chunk_size, overlap, test_chunks_file) == 0);

    // Run fasta merge to merge back a combined sequence file
    CuAssertTrue(testCase, st_system("faffy merge --logLevel DEBUG -i %s -o %s", test_chunks_file, test_dechunked_fa_file) == 0);

    // Check the sequences are equal
    FILE *fh = fopen(test_fa_file, "r"); stHash *seqs = fastaReadToMap(fh); fclose(fh);
    fh = fopen(test_dechunked_fa_file, "r"); stHash *seqs2 = fastaReadToMap(fh); fclose(fh);
    CuAssertIntEquals(testCase, stHash_size(seqs), stHash_size(seqs2));
    stHashIterator *it = stHash_getIterator(seqs);
    char *header;
    while((header = stHash_getNext(it)) != NULL) {
        CuAssertTrue(testCase, stHash_search(seqs, header) != NULL);
        CuAssertTrue(testCase, stHash_search(seqs2, header) != NULL);
        CuAssertIntEquals(testCase, strlen(stHash_search(seqs, header)), strlen(stHash_search(seqs2, header)));
        CuAssertTrue(testCase, strcmp(stHash_search(seqs, header), stHash_search(seqs2, header)) == 0);
    }
    stHash_destructIterator(it);
    stHash_destruct(seqs);
    stHash_destruct(seqs2);

    // Cleanup
    CuAssertTrue(testCase, st_system("rm -rf %s %s %s %s", test_fa_file, test_fasta_chunks_dir, test_chunks_file, test_dechunked_fa_file) == 0);
}

CuSuite* addFastaChunkAndMergeTestSuite(void) {
    CuSuite* suite = CuSuiteNew();
    SUITE_ADD_TEST(suite, test_fasta_chunk_and_merge);
    return suite;
}
