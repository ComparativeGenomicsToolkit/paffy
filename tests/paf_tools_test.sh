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

# Run paffy to_bed basic
echo "paffy to_bed basic output"
paffy to_bed -i ${working_dir}/output.paf -o ${working_dir}/output.bed
[ -s ${working_dir}/output.bed ]

# Run paffy to_bed -b (binary: 4th column must be 0 or 1)
echo "paffy to_bed binary output"
paffy to_bed -i ${working_dir}/output.paf -b -o ${working_dir}/output_binary.bed
[ -s ${working_dir}/output_binary.bed ]
awk '{if ($4 > 1) exit 1}' ${working_dir}/output_binary.bed

# Run paffy to_bed -e (exclude unaligned: all output rows must have count > 0)
echo "paffy to_bed exclude unaligned"
paffy to_bed -i ${working_dir}/output.paf -e -o ${working_dir}/output_no_unaligned.bed
[ -s ${working_dir}/output_no_unaligned.bed ]
awk '{if ($4 == 0) exit 1}' ${working_dir}/output_no_unaligned.bed

# Run paffy to_bed -f (exclude aligned: all output rows must have count == 0)
echo "paffy to_bed exclude aligned"
paffy to_bed -i ${working_dir}/output.paf -f -o ${working_dir}/output_unaligned_only.bed
awk '{if ($4 != 0) exit 1}' ${working_dir}/output_unaligned_only.bed

# Run paffy to_bed -n (include inverted: flips alignments so target seqs are also covered)
echo "paffy to_bed include inverted"
lines_non_inv=$(paffy to_bed -i ${working_dir}/output.paf -e | wc -l)
lines_inv=$(paffy to_bed -i ${working_dir}/output.paf -e -n | wc -l)
[ "${lines_inv}" -ge "${lines_non_inv}" ]

# Run paffy filter -u (min identity after encoding mismatches)
echo "paffy filter by min identity"
paffy add_mismatches -i ${working_dir}/output.paf ${working_dir}/*.fa \
  | paffy filter -u 0.7 > /dev/null

# Run paffy filter -v (min identity with gaps after encoding mismatches)
echo "paffy filter by min identity with gaps"
paffy add_mismatches -i ${working_dir}/output.paf ${working_dir}/*.fa \
  | paffy filter -v 0.7 > /dev/null

# Run paffy filter -w (max tile level after tile)
echo "paffy filter by max tile level"
paffy tile -i ${working_dir}/output.paf | paffy filter -w 1 > /dev/null

# Run paffy trim -f (fixed trim: constant fraction from each end)
# Fixed trim reduces total aligned bases and slightly lowers identity, so use a relaxed -u threshold
echo "paffy trim fixed trim"
paffy trim -f -t 0.1 -i ${working_dir}/output.paf | paffy view ${working_dir}/*.fa -s -t -u 0.73 -v 450000

# Run paffy dedupe -a (check inverse: also deduplicate inverted alignments)
echo "paffy dedupe check inverse"
paffy invert -i ${working_dir}/output.paf > ${working_dir}/output_inv.paf
cat ${working_dir}/output.paf ${working_dir}/output_inv.paf \
  | paffy dedupe -a | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000
