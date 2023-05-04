#!/usr/bin/env bash

# Script to test pafs

# exit when any command fails
set -e

# log the commands its running
set -x

# Make a working directory
#working_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
working_dir=./temp_tools
mkdir -p ${working_dir}

# Make sure we cleanup the temp dir
trap "rm -rf ${working_dir}" EXIT

# Get the sequences
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6 -O ${working_dir}/simCow.chr6.fa
wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6 -O ${working_dir}/simDog.chr6.fa

# Run lastz
lastz ${working_dir}/*.fa --format=paf > ${working_dir}/output.paf

# Run paffy view
echo "minimum local alignment identity"
paffy view -i ${working_dir}/output.paf ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with invert
echo "paffy invert minimum local alignment identity"
paffy invert -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with chain
echo "paffy chain minimum local alignment identity"
paffy chain -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with shatter
echo "paffy shatter minimum local alignment identity (will be low as equal to worst run of matches)"
paffy shatter -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with tile
echo "paffy tile minimum local alignment identity"
paffy tile -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy add_mismatches
echo "adding mismatches"
paffy add_mismatches -i ${working_dir}/output.paf ${working_dir}/*.fa | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy add_mismatches
echo "adding and then remove mismatches"
paffy add_mismatches -i ${working_dir}/output.paf ${working_dir}/*.fa | paffy add_mismatches -a | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with trim (identity may be higher as we trim the tails)
echo "paffy trim minimum local alignment identity"
paffy add_mismatches -i ${working_dir}/output.paf ${working_dir}/*.fa | paffy trim -r 0.05 | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 446000

# Run paffy view with trim (identity may be higher as we trim the tails)
echo "paffy trim minimum local alignment identity, ignoring mismatches"
paffy trim -r 0.95 -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy view with filter
echo "paffy filter removing alignments with score < than 10000"
paffy filter -i ${working_dir}/output.paf -t 10000 | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 515000

# Run paffy view with filter (inverted)
echo "paffy filter removing alignments with score >= than 10000"
paffy filter -i ${working_dir}/output.paf -t 10000 -x | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 15000
