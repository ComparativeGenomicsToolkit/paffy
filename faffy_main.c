/*
 * faffy: The faffy command line interface, whose subcommands run different tools
 *
 *  Released under the MIT license, see LICENSE.txt
*/

#include "paf.h"
#include <getopt.h>
#include <time.h>

extern int fasta_chunk_main(int argc, char *argv[]);
extern int fasta_extract_main(int argc, char *argv[]);
extern int fasta_merge_main(int argc, char *argv[]);

void usage(void) {
    fprintf(stderr, "faffy: little toolkit for working with FASTA files\n\n");
    fprintf(stderr, "usage: faffy <command> [options]\n\n");
    fprintf(stderr, "available commands:\n");
    fprintf(stderr, "    chunk                  Break a large fasta file into smaller files for parallel processing\n");
    fprintf(stderr, "    merge                  Merge together the chunks created by chunk, potentially resolving overlaps\n");
    fprintf(stderr, "    extract                Extract subsequences of the fasta file\n");
    fprintf(stderr, "\n");
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        usage();
        return 0;
    }

    if (strcmp(argv[1], "chunk") == 0) {
        return fasta_chunk_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "merge") == 0) {
        return fasta_merge_main(argc - 1, argv + 1);
    } else if (strcmp(argv[1], "extract") == 0) {
        return fasta_extract_main(argc - 1, argv + 1);
    } else {
        fprintf(stderr, "%s is not a valid faffy command\n", argv[1]);
        usage();
        return 1;
    }
    return 1;
}
