# A small C/CLI library for manipulating PAF files.

This is a C and CLI library for manipulating/reading/writing 
[Pairwise Alignment Format (PAF)](https://github.com/lh3/miniasm/blob/master/PAF.md).

# Installing Paffy CLI/C Library

Do build this repo clone the repo as follows and then make:

    git clone https://github.com/ComparativeGenomicsToolkit/paffy.git --recursive
    cd paf && make

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
    tile           Give alignments levels, from lowest (best) to highest (worse) by greedily picking
                   the best alignment at each location
    dedupe         Remove duplicate alignments from a file based on exact query/target coordinates
    dechunk        Manipulate coordinates to allow aggregation of PAFs computed over subsequences
    upconvert      Converts the coordinates of paf alignments to refer to extracted subsequences
```

For example, to pretty print a PAF alignment:

    paffy view -i PAF_FILE [FASTA_FILES] -a

# Using C Library

There is also a simple C library for working with taf/maf files. See paf.h in the
inc directory.

