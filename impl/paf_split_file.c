/*
 * paffy split_file: Split PAF file into per-contig output files
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include <getopt.h>
#include <time.h>

static void usage(void) {
     fprintf(stderr, "paffy split_file [options], version 0.1\n");
     fprintf(stderr, "Split PAF file into separate output files by target (default) or query contig name\n");
     fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
     fprintf(stderr, "-p --prefix : Output file prefix (may include directory path). Default: split_\n");
     fprintf(stderr, "-q --query : Split by query contig name instead of target contig name\n");
     fprintf(stderr, "-m --minLength : Contigs with sequence length < m are grouped into combined files\n"
                     "                 (<prefix>small_0.paf, <prefix>small_1.paf, ...) such that the total\n"
                     "                 contig length in each file does not exceed m. All alignments for a\n"
                     "                 given contig go in exactly one file. Default: 0 (disabled)\n");
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
     bool split_by_query = 0;
     int64_t minLength = 0;

     ///////////////////////////////////////////////////////////////////////////
     // Parse the inputs
     ///////////////////////////////////////////////////////////////////////////

     while (1) {
         static struct option long_options[] = { { "logLevel", required_argument, 0, 'l' },
                                                 { "inputFile", required_argument, 0, 'i' },
                                                 { "prefix", required_argument, 0, 'p' },
                                                 { "query", no_argument, 0, 'q' },
                                                 { "minLength", required_argument, 0, 'm' },
                                                 { "help", no_argument, 0, 'h' },
                                                 { 0, 0, 0, 0 } };

         int option_index = 0;
         int64_t key = getopt_long(argc, argv, "l:i:p:qm:h", long_options, &option_index);
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
             case 'q':
                 split_by_query = 1;
                 break;
             case 'm':
                 minLength = atol(optarg);
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
     st_logInfo("Split by : %s\n", split_by_query ? "query" : "target");
     st_logInfo("Min contig length : %" PRIi64 "\n", minLength);

     //////////////////////////////////////////////
     // Split the paf file
     //////////////////////////////////////////////

     FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
     stHash *contig_to_file = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);

     // For small contigs: map contig_name -> FILE* so all alignments for a contig go to the same file
     stHash *small_contig_to_file = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, free, NULL);
     stList *small_files = stList_construct(); // list of FILE* for closing
     FILE *current_small_file = NULL;
     int64_t current_small_file_length = 0;
     int64_t small_file_index = 0;

     Paf *paf;
     int64_t total_records = 0;
     while((paf = paf_read(input, 0)) != NULL) {
         char *contig_name = split_by_query ? paf->query_name : paf->target_name;
         int64_t contig_length = split_by_query ? paf->query_length : paf->target_length;
         FILE *output;
         if (minLength > 0 && contig_length < minLength) {
             // Check if this small contig already has an assigned file
             output = stHash_search(small_contig_to_file, contig_name);
             if (output == NULL) {
                 // New small contig - check if it fits in the current small file
                 if (current_small_file == NULL || current_small_file_length + contig_length > minLength) {
                     // Start a new small file
                     char *filename = stString_print("%ssmall_%" PRIi64 ".paf", prefix, small_file_index++);
                     current_small_file = fopen(filename, "w");
                     if (current_small_file == NULL) {
                         st_errAbort("Could not open output file: %s\n", filename);
                     }
                     st_logInfo("Opened small contigs output file: %s\n", filename);
                     free(filename);
                     stList_append(small_files, current_small_file);
                     current_small_file_length = 0;
                 }
                 current_small_file_length += contig_length;
                 stHash_insert(small_contig_to_file, stString_copy(contig_name), current_small_file);
                 output = current_small_file;
             }
         } else {
             output = get_output_file(contig_to_file, contig_name, prefix);
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

     // Close all per-contig output files
     stHashIterator *it = stHash_getIterator(contig_to_file);
     char *key;
     while ((key = stHash_getNext(it)) != NULL) {
         FILE *fh = stHash_search(contig_to_file, key);
         fclose(fh);
     }
     stHash_destructIterator(it);
     stHash_destruct(contig_to_file);

     // Close all small contig output files
     for (int64_t i = 0; i < stList_length(small_files); i++) {
         fclose(stList_get(small_files, i));
     }
     stList_destruct(small_files);
     stHash_destruct(small_contig_to_file);

     st_logInfo("Paffy split_file is done! Split %" PRIi64 " records, %" PRIi64 " seconds have elapsed\n",
                total_records, time(NULL) - startTime);

     return 0;
 }
