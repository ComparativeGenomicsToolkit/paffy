/*
 * paffy filter: Filter based on alignment stats
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include <getopt.h>
#include <time.h>

static void usage(void) {
     fprintf(stderr, "paffy filter [options], version 0.1\n");
     fprintf(stderr, "Filter pafs based on alignment stats\n");
     fprintf(stderr, "-i --inputFile : Input paf file. If not specified reads from stdin\n");
     fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
     fprintf(stderr, "-s --minChainScore : Filter alignments with a chain score less than this\n");
     fprintf(stderr, "-t --minAlignmentScore : Filter alignments with an alignment score less than this\n");
     fprintf(stderr, "-u --minIdentity : Filter alignments with an identity less than this, exclude indels\n");
     fprintf(stderr, "-v --minIdentityWithGaps : Filter alignments with an identity less than this, including indels\n");
     fprintf(stderr, "-w --maxTileLevel : Filter alignments with a tile level greater than this\n");
     fprintf(stderr, "-x --invert : Only output alignments that don't pass filters\n");
     fprintf(stderr, "-l --logLevel : Set the log level\n");
     fprintf(stderr, "-h --help : Print this help message\n");
 }

 int paffy_filter_main(int argc, char *argv[]) {
     time_t startTime = time(NULL);
     int64_t min_chain_score = -1;
     int64_t min_alignment_score = -1;
     double min_identity = -1.0;
     double min_identity_with_gaps = -1.0;
     int64_t max_tile_level = -1;
     bool invert = 0;

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
                                                 { "minChainScore", required_argument, 0, 's' },
                                                 { "minAlignmentScore", required_argument, 0, 't' },
                                                 { "minIdentity", required_argument, 0, 'u' },
                                                 { "minIdentityWithGaps", required_argument, 0, 'v' },
                                                 { "maxTileLevel", required_argument, 0, 'w' },
                                                 { "invert", no_argument, 0, 'x' },
                                                 { "help", no_argument, 0, 'h' },
                                                 { 0, 0, 0, 0 } };

         int option_index = 0;
         int64_t key = getopt_long(argc, argv, "l:i:o:s:t:u:v:w:xh", long_options, &option_index);
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
             case 's':
                 min_chain_score = atoi(optarg);
                 break;
             case 't':
                 min_alignment_score = atoi(optarg);
                 break;
             case 'u':
                 min_identity = atof(optarg);;
                 break;
             case 'v':
                 min_identity_with_gaps = atof(optarg);;
                 break;
             case 'w':
                 max_tile_level = atoi(optarg);
                 break;
             case 'x':
                 invert = 1;
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
     st_logInfo("Filtering paf with min chain score:%" PRIi64 " min alignment score:%" PRIi64
                " min identity:%f min identity with gaps:%f max tile level:%" PRIi64 " invert:%s\n", min_chain_score,
                min_alignment_score, min_identity, min_identity_with_gaps, max_tile_level, invert ? "True" : "False");

     //////////////////////////////////////////////
     // Filter the paf
     //////////////////////////////////////////////

     FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
     FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

     Paf *paf;
     while((paf = paf_read(input)) != NULL) {
         // Calculate identity stats
         int64_t matches=0, mismatches=0, query_inserts=0, query_deletes=0,
                 query_insert_bases=0, query_delete_bases=0;
         paf_stats_calc(paf, &matches, &mismatches, &query_inserts,
                        &query_deletes, &query_insert_bases, &query_delete_bases, 0);
         double identity = (float)matches / (matches + mismatches);
         double identity_with_gaps = (float)matches / (matches + mismatches + query_insert_bases + query_delete_bases);
         if(paf->score >= min_alignment_score && paf->chain_score >= min_chain_score &&
            (max_tile_level == -1 || paf->tile_level <= max_tile_level) && identity >= min_identity &&
            identity_with_gaps >= min_identity_with_gaps) {
             if(invert) {
                 if(st_getLogLevel() == debug) {
                     st_logDebug("Filtering alignment with matches:%" PRIi64 ", identity: %f (%f with gaps), score: %" PRIi64
                     ", chain-score:%" PRIi64 "\n", matches, identity, identity_with_gaps, paf->score, paf->chain_score);
                     paf_write(paf, stderr);
                 }
             }
             else {
                 paf_write(paf, output);
             }
         }
         else {
             if(invert) {
                 paf_write(paf, output);
             }
             else if(st_getLogLevel() == debug) {
                st_logDebug("Filtering alignment with matches:%" PRIi64 ", identity: %f (%f with gaps), score: %" PRIi64
                 ", chain-score:%" PRIi64 "\n", matches, identity, identity_with_gaps, paf->score, paf->chain_score);
                paf_write(paf, stderr);
            }
         }
         paf_destruct(paf);
     }

     //////////////////////////////////////////////
     // Cleanup
     //////////////////////////////////////////////

     if(inputFile != NULL) {
         fclose(input);
     }
     if(outputFile != NULL) {
         fclose(output);
     }

     st_logInfo("Paffy filter is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

     //while(1);
     //assert(0);

     return 0;
 }
