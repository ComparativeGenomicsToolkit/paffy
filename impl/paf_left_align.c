/*
 * paffy left_align: Move each gap of a paf's alignment as far left on the target as the sequences allow
 *
 *  Released under the MIT license, see LICENSE.txt
 *
 * An indel in a repeat can be written anywhere along the repeat with the same score, and aligners choose
 * between those places inconsistently (minigraph, for one, pins them to seed anchors, and puts them at
 * mirrored places on the two strands). Aligning every gap to the left of the target's forward strand
 * gives alignments to the same target one representation of each indel.
 *
 * Overview:
 * (1) Load query and target sequences
 * (2) For each input PAF record left-align the gaps
*/

#include "paf.h"
#include <getopt.h>
#include <time.h>
#include "bioioC.h"

static void usage(void) {
    fprintf(stderr, "paffy left_align [fasta_files]xN [options], version 0.1\n");
    fprintf(stderr, "Move each gap as far left on the target as the sequences allow, without changing the alignment's score\n");
    fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
    fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
    fprintf(stderr, "-l --logLevel : Set the log level\n");
    fprintf(stderr, "-h --help : Print this help message\n");
}

int paffy_left_align_main(int argc, char *argv[]) {
    time_t startTime = time(NULL);

    /*
     * Arguments/options
     */
    char *logLevelString = NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;

    ///////////////////////////////////////////////////////////////////////////
    // Parse the inputs
    ///////////////////////////////////////////////////////////////////////////

    while (1) {
        static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                { "inputFile", required_argument, 0, 'i' },
                                                { "outputFile", required_argument, 0, 'o' },
                                                { "help", no_argument, 0, 'h' },
                                                { 0, 0, 0, 0 } };

        int option_index = 0;
        int64_t key = getopt_long(argc, argv, "l:i:o:h", long_options, &option_index);
        if (key == -1) {
            break;
        }

        switch (key) {
            case 'l':
                logLevelString = optarg;
                break;
            case 'i':
                inputFile = optarg;
                break;
            case 'o':
                outputFile = optarg;
                break;
            case 'h':
                usage();
                return 0;
            default:
                usage();
                return 1;
        }
    }

    //////////////////////////////////////////////
    //Log the inputs
    //////////////////////////////////////////////

    st_setLogLevelFromString(logLevelString);
    st_logInfo("Input file string : %s\n", inputFile);
    st_logInfo("Output file string : %s\n", outputFile);

    //////////////////////////////////////////////
    // Parse the sequences
    //////////////////////////////////////////////

    stHash *sequences = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, free);
    while(optind < argc) {
        char *seq_file = argv[optind++];
        st_logInfo("Parsing sequence file : %s\n", seq_file);
        FILE *seq_file_handle = fopen(seq_file, "r");
        if(seq_file_handle == NULL) {
            st_errAbort("Could not open sequence file: %s\n", seq_file);
        }
        fastaReadToFunction(seq_file_handle, sequences, fastaRead_readToMapFunction);
        fclose(seq_file_handle);
    }
    st_logInfo("Read %i sequences from sequence files\n", (int)stHash_size(sequences));

    //////////////////////////////////////////////
    // Left-align the paf records
    //////////////////////////////////////////////

    FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
    FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");
    if(input == NULL || output == NULL) {
        st_errAbort("Could not open %s\n", input == NULL ? inputFile : outputFile);
    }

    int64_t records = 0, gaps_moved = 0;
    Paf *paf;
    int64_t paf_buffer_length = 100;
    char *paf_buffer = st_malloc(sizeof(char) * paf_buffer_length);
    while((paf = paf_read_with_buffer(input, 1, &paf_buffer, &paf_buffer_length)) != NULL) {
        if(paf->cigar != NULL) {
            char *query_seq = stHash_search(sequences, paf->query_name);
            if(query_seq == NULL) {
                fprintf(stderr, "No query sequence named: %s found\n", paf->query_name);
                exit(1);
            }
            char *target_seq = stHash_search(sequences, paf->target_name);
            if(target_seq == NULL) {
                fprintf(stderr, "No target sequence named: %s found\n", paf->target_name);
                exit(1);
            }
            gaps_moved += paf_left_align(paf, query_seq, target_seq);
        }

        // Check all is good
        paf_check(paf);

        paf_write_with_buffer(paf, output, &paf_buffer, &paf_buffer_length);

        paf_destruct(paf);
        records++;
    }
    free(paf_buffer);

    //////////////////////////////////////////////
    // Cleanup
    //////////////////////////////////////////////

    if(inputFile != NULL) {
        fclose(input);
    }
    if(outputFile != NULL) {
        st_fclose(output, outputFile);
    }
    stHash_destruct(sequences);

    st_logInfo("Moved %" PRIi64 " gaps in %" PRIi64 " records\n", gaps_moved, records);
    st_logInfo("Paffy left_align is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

    return 0;
}
