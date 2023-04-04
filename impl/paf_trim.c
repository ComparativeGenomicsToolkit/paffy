/*
 * paf_trim: Trim bases from the prefix and suffix matches of a paf alignment
 *   (1) Load local alignments file (PAF)
 *   (2) Trim prefix and suffix alignment, either a constant amount or according to an identity threshold.
 *   (3) Output local alignments file (PAF)
 *
 * Released under the MIT license, see LICENSE.txt
 */

#include "paf.h"
#include <getopt.h>
#include <time.h>

static float trim_end_fraction = 1.0;
static float trim_by_identity = 1;
static float trim_by_identity_fraction = 0.3;

static void usage(void) {
     fprintf(stderr, "paf_trim [options], version 0.1\n");
     fprintf(stderr, "Trims the ends of a PAF file\n");
     fprintf(stderr, "-i --inputFile : Input paf file to invert. If not specified reads from stdin\n");
     fprintf(stderr, "-o --outputFile : Output paf file. If not specified outputs to stdout\n");
     fprintf(stderr, "-r --trimIdentity : Trim tails with\n"
                     "alignment identity lower than this fraction of the overall alignment identity (from 0 to 1,\n"
                     "by default: %f). If mismatches are not encoded in the cigar then identity is fraction of aligned\n"
                     "bases, if mismatches in are encoded identity is fraction of aligned and matched bases.\n", trim_by_identity_fraction);
     fprintf(stderr, "-t --trimFraction : Fraction (from 0 to 1) of aligned bases to trim from each end of the \n"
                     "alignment (default:%f). If not --fixedTrim (see below) trimFraction is the \n"
                     "max fraction of the alignment to trim in each tail\n", trim_end_fraction);
     fprintf(stderr, "-f --fixedTrim : Trim a constant amount from each tail instead of trimming by identity. Amount\n"
                     " to trim is determined by --trimFraction\n");
     fprintf(stderr, "-l --logLevel : Set the log level\n");
     fprintf(stderr, "-h --help : Print this help message\n");
 }

 int paf_trim_main(int argc, char *argv[]) {
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
                                                 { "trimFraction", required_argument, 0, 't' },
                                                 { "trimIdentity", required_argument, 0, 'r' },
                                                 { "fixedTrim", no_argument, 0, 'f' },
                                                 { "help", no_argument, 0, 'h' },
                                                 { 0, 0, 0, 0 } };

         int option_index = 0;
         int64_t key = getopt_long(argc, argv, "l:i:o:ht:r:f", long_options, &option_index);
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
             case 't':
                 trim_end_fraction = atof(optarg);
                 break;
             case 'r':
                 trim_by_identity_fraction = atof(optarg);
                 break;
             case 'f':
                 trim_by_identity = 0;
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
     st_logInfo("Trim fraction using : %f\n", trim_end_fraction);
     st_logInfo("Trim by identity fraction : %f\n", trim_by_identity_fraction);

     //////////////////////////////////////////////
     // Invert the paf
     //////////////////////////////////////////////

     FILE *input = inputFile == NULL ? stdin : fopen(inputFile, "r");
     FILE *output = outputFile == NULL ? stdout : fopen(outputFile, "w");

     Paf *paf;
     while((paf = paf_read(input)) != NULL) {
         if(trim_by_identity) {
             paf_trim_unreliable_tails(paf, trim_by_identity_fraction, trim_end_fraction);
         }
         else {
             paf_trim_end_fraction(paf, trim_end_fraction); // the invert routine
         }

         paf_check(paf);
         paf_write(paf, output);
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

     st_logInfo("Paf invert is done!, %" PRIi64 " seconds have elapsed\n", time(NULL) - startTime);

     //while(1);
     //assert(0);

     return 0;
 }
