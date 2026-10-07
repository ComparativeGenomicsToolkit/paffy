# A small C/CLI library for manipulating PAF files.

This is a C and CLI library for manipulating/reading/writing 
[Pairwise Alignment Format (PAF)](https://github.com/lh3/miniasm/blob/master/PAF.md).
It is used by [Cactus](https://github.com/ComparativeGenomicsToolkit/cactus) in creating the 
local, pairwise alignments between input sequences. It also contains some utilities
for manipulating [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/#:~:text=In%20FASTA%20format%20the%20line,should%20not%20contain%20any%20spaces.) 
files.

# Installing Paffy CLI/C Library

Do build this repo clone the repo as follows and then make:

    git clone https://github.com/ComparativeGenomicsToolkit/paffy.git --recursive
    cd paffy && make

To test the installation, after adding the paffy/bin directory to your path, do:

    make test

This will run a small number of tests. 

# Paffy Utilities

All Paffy utilities are run using `paffy <command>`, where the available commands are:

```
    view           Pretty print and extract stats about PAF alignments
    chain          Chain together paf alignments, useful for getting synteny info
    add_mismatches Add mismatch information to the PAF cigars (i.e. convert M to =/X format)
    trim           Slice of lower identity tail alignments
    to_bed         Build an alignment coverage map of a chosen sequence in BED format
    shatter        Break the PAFs into gapless subalignments
    invert         Switch query and target
    left_align     Move each gap along its repeat (score unchanged), so equivalent indels are written the same
                   way in every alignment: left on the target, or with --canonical at whichever end the repeat
                   reads smaller than its reverse complement, which agrees across any pairs and strands
    tile           Give alignments levels, from lowest (best) to highest (worse) by greedily picking
                   the best alignment at each location
    dedupe         Remove duplicate alignments from a file based on exact query/target coordinates
    dechunk        Manipulate coordinates to allow aggregation of PAFs computed over subsequences
    upconvert      Converts the coordinates of paf alignments to refer to extracted subsequences
    split_file     Split a PAF file into separate output files by target contig name. Optionally
                   group small contigs (below a given target length threshold) into size-bounded files
    unanchor       Remove the minigraph (hub) anchors of every genome at tandem-repeat loci where the genomes'
                   anchors disagree, so that cactus realigns those loci itself (see below)
```

`paffy unanchor` takes one chromosome's PAF of genomes against minigraph nodes (cactus-align's, after
filter_paf), the node FASTA and the minigraph rGFA:

    paffy unanchor -i chr.paf -n chr.gfa.fa.gz -g chr.gfa.gz -r CHM13 -o chr.cut.paf -L loci.bed -b hub.bed -v veto.bed -q query.bed

It tandem-scans the reference genome's path through the graph (periods 2-100, satellite arrays excluded), walks
each other contig's anchors along that path, and picks the repeat loci where at least `--minHaps` contigs both
insert and delete bases between the anchors on either side (they disagree on the repeat's copy number or
structure). A locus is not picked when something would make clearing it unsafe: it spans more than half of
cactus's banding limit, the reference does not anchor itself there, a contig's anchors cross its edge, flip
strand or step back inside it, an alt node through it is used elsewhere, some genome's sequence between the
flanking anchors is too long, an alt allele through it holds non-repeat sequence that is not from the locus (an
SV), or it has neither an alt node nor a node boundary. Every matched column on the picked loci's nodes is
removed from every record, the reference's included, and the query bases that lost all their anchors are written
to `--queryBed`. The `--logFile` summarises every step; the other BEDs list the cleared loci, the cleared node
intervals and the loci not cleared with their reasons. The output does not depend on `--threads`.

In addition the FASTA utilities are run using the `faffy <command>`, where the available commands are:
```
    chunk           Break a fasta file into smaller files for parallel computation
    merge           Merge together fasta files, dealing with overlap (the inverse of chunk)
    extract         Extract subsequences from a fasta file 
```

For example, to pretty print a PAF alignment:

    paffy view -i PAF_FILE [FASTA_FILES] -a

To understand the chunking/dechunking code for fastas and pafs see the
tests/faf_paf_chunking_test.sh.

# Using C Library

There is also a simple C library for working with paf files. See paf.h in the
inc directory.

# Cactus Local Alignment Pipeline

To understand the steps in the pairwise alignment pipeline used by Cactus
see tests/paf_pipeline_test.sh. The objective is to find chains
of high-quality local alignments used to anchor the multiple sequence alignment
process.


# Todo

Create a clean script for reproducing Cactus local alignment pipeline

Improve docs
