/*
 * paffy: The paffy command line interface, whose subcommands run different tools
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "paf.h"
#include <getopt.h>
#include <time.h>

extern int paffy_add_mismatches_main(int argc, char *argv[]);
extern int paffy_chain_main(int argc, char *argv[]);
extern int paffy_dechunk_main(int argc, char *argv[]);
extern int paffy_dedupe_main(int argc, char *argv[]);
extern int paffy_invert_main(int argc, char *argv[]);
extern int paffy_shatter_main(int argc, char *argv[]);
extern int paffy_tile_main(int argc, char *argv[]);
extern int paffy_to_bed_main(int argc, char *argv[]);
extern int paffy_trim_main(int argc, char *argv[]);
extern int paffy_upconvert_main(int argc, char *argv[]);
extern int paffy_view_main(int argc, char *argv[]);
extern int paffy_filter_main(int argc, char *argv[]);

void usage(void) {
    fprintf(stderr, "paffy: toolkit for working with PAF files\n\n");
    fprintf(stderr, "usage: paffy <command> [options]\n\n");
    fprintf(stderr, "available commands:\n");
    fprintf(stderr, "    add_mismatches           Replace Ms with =/Xs in PAF cigar string\n");
    fprintf(stderr, "    chain                    Chain together PAF alignments\n");
    fprintf(stderr, "    dechunk                  Manipulate coordinates to allow aggregation of PAFs computed over subsequences\n");
    fprintf(stderr, "    dedupe                   Remove duplicate alignments from a file based on exact query/target coordinates\n");
    fprintf(stderr, "    filter                   Filter alignments based upon alignment stats\n");
    fprintf(stderr, "    invert                   Switch query and target coordinates\n");
    fprintf(stderr, "    shatter                  Break PAFs into sequence of gapless PAF alignments\n");
    fprintf(stderr, "    tile                     Give alignments levels, from lowest (best) to highest (worse) by greedily picking\n"
                    "                             the best alignment at each location\n");
    fprintf(stderr, "    to_bed                   Build an alignment coverage map of a chosen sequence in BED format\n");
    fprintf(stderr, "    trim                     Slice of lower identity tail alignments\n");
    fprintf(stderr, "    upconvert                Converts the coordinates of paf alignments to refer to extracted subsequences\n");
    fprintf(stderr, "    view                     Pretty print and extract stats about PAF alignments\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage();
        return 0;
    }

    if (strcmp(argv[1], "add_mismatches") == 0) {
        return paffy_add_mismatches_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "chain") == 0) {
        return paffy_chain_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "dechunk") == 0) {
        return paffy_dechunk_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "dedupe") == 0) {
        return paffy_dedupe_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "invert") == 0) {
        return paffy_invert_main(argc - 1, argv + 1);
    }else if (strcmp(argv[1], "filter") == 0) {
        return paffy_filter_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "shatter") == 0) {
        return paffy_shatter_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "tile") == 0) {
        return paffy_tile_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "to_bed") == 0) {
        return paffy_to_bed_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "trim") == 0) {
        return paffy_trim_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "upconvert") == 0) {
        return paffy_upconvert_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "view") == 0) {
        return paffy_view_main(argc - 1, argv + 1);
    } else {
        fprintf(stderr, "%s is not a valid paffy command\n", argv[1]);
        usage();
        return 1;
    }
    return 1;
}
