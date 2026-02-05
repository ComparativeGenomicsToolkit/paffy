/*
 * paffy split_file: Split PAF file into per-target-contig output files
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include <getopt.h>
#include <time.h>

static void usage(void) {
     fprintf(stderr, "paffy split_file [options], version 0.1\n");
     fprintf(stderr, "Split PAF file into separate output files by target contig name\n");
     fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
     fprintf(stderr, "-p --prefix : Output file prefix (may include directory path). Default: split_\n");
     fprintf(stderr, "-m --minTargetLength : Contigs with target_length < m go to <prefix>small.paf. Default: 0 (disabled)\n");
     fprintf(stderr, "-l --logLevel : Set the log level\n");
     fprintf(stderr, "-h --help : Print this help message\n");
 }

/*
 * Sanitize a target name for use in a filename by replacing '/' with '_'
 */
static char *sanitize_filename(const char *name) {
    char *sanitized = stString_copy(name);
    for (int64_t i = 0; sanitized[i] != '\0'; i++) {
        if (sanitized[i] == '/') {
            sanitized[i] = '_';
        }
    }
    return sanitized;
}

/*
 * Get or create an output file handle for the given target name
 */
static FILE *get_output_file(stHash *target_to_file, const char *target_name, const char *prefix) {
    FILE *fh = stHash_search(target_to_file, (void *)target_name);
    if (fh == NULL) {
        char *sanitized = sanitize_filename(target_name);
        char *filename = stString_print("%s%s.paf", prefix, sanitized);
        fh = fopen(filename, "w");
        if (fh == NULL) {
            st_errAbort("Could not open output file: %s\n", filename);
        }
        st_logInfo("Opened output file: %s\n", filename);
        stHash_insert(target_to_file, stString_copy(target_name), fh);
        free(sanitized);
        free(filename);
    }
    return fh;
}

 int paffy_split_file_main(int argc, char *argv[]) {
     time_t startTime = time(NULL);

     /*
      * Arguments/options
      */
     char *logLevelString = NULL;
     char *inputFile = NULL;
     char *prefix = "split_";
     int64_t minTargetLength = 0;

     ///////////////////////////////////////////////////////////////////////////
     // Parse the inputs
     ///////////////////////////////////////////////////////////////////////////

     while (1) {
         static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                 { "inputFile", required_argument, 0, 'i' },
                                                 { "prefix", required_argument, 0, 'p' },
                                                 { "minTargetLength", required_argument, 0, 'm' },
                                                 { "help", no_argument, 0, 'h' },
                                                 { 0, 0, 0, 0 } };

         int option_index = 0;
         int64_t key = getopt_long(argc, argv, "l:i:p:m:h", long_options, &option_index);
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
             case 'p':
                 prefix = optarg;
                 break;
             case 'm':
                 minTargetLength = atol(optarg);
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
     st_logInfo("Output prefix : %s\n", prefix);
     st_logInfo("Min target length : %" PRIi64 "\n", minTargetLength);

     //////////////////////////////////////////////
     // Split the paf file
     //////////////////////////////////////////////

     FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
     stHash *target_to_file = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
     FILE *small_file = NULL;

     Paf *paf;
     int64_t total_records = 0;
     while((paf = paf_read(input, 0)) != NULL) {
         FILE *output;
         if (minTargetLength > 0 && paf->target_length < minTargetLength) {
             // Lazily open the small contigs file
             if (small_file == NULL) {
                 char *filename = stString_print("%ssmall.paf", prefix);
                 small_file = fopen(filename, "w");
                 if (small_file == NULL) {
                     st_errAbort("Could not open output file: %s\n", filename);
                 }
                 st_logInfo("Opened small contigs output file: %s\n", filename);
                 free(filename);
             }
             output = small_file;
         } else {
             output = get_output_file(target_to_file, paf->target_name, prefix);
         }
         paf_write(paf, output);
         total_records++;
         paf_destruct(paf);
     }

     //////////////////////////////////////////////
     // Cleanup
     //////////////////////////////////////////////

     if(inputFile != NULL) {
         fclose(input);
     }

     // Close all per-target output files
     stHashIterator *it = stHash_getIterator(target_to_file);
     char *key;
     while ((key = stHash_getNext(it)) != NULL) {
         FILE *fh = stHash_search(target_to_file, key);
         fclose(fh);
     }
     stHash_destructIterator(it);
     stHash_destruct(target_to_file);

     if (small_file != NULL) {
         fclose(small_file);
     }

     st_logInfo("Paffy split_file is done! Split %" PRIi64 " records, %" PRIi64 " seconds have elapsed\n",
                total_records, time(NULL) - startTime);

     return 0;
 }
