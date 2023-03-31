#!/usr/bin/env bash

# Script to test pafs

# exit when any command fails
set -e

# log the commands its running
set -x

# Make a working directory
#working_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
working_dir=./temp_chains
mkdir ${working_dir}

# Make sure we cleanup the temp dir
#trap "rm -rf ${working_dir}" EXIT

# Get the sequences
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanPan1_XY_1_5000000.fa -O ${working_dir}/mPanPan1_XY_1_5000000.fa
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanTro3_XY_1_5000000.fa -O ${working_dir}/mPanTro3_XY_1_5000000.fa

# Run lastz
lastz ${working_dir}/mPanPan1_XY_1_5000000.fa[multiple][nameparse=darkspace] ${working_dir}/mPanTro3_XY_1_5000000.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz.paf

# Inverting so for every query:target alignment we have the mirror target:query alignment, this is used to
# make alignments symmetric for chaining
echo "Inverting"
paf_invert -i ${working_dir}/lastz.paf > ${working_dir}/inverted.paf

# Catting
echo "Catting forward and inverted pafs"
cat ${working_dir}/lastz.paf ${working_dir}/inverted.paf > ${working_dir}/combined.paf

# Run paf_add_mismatches
echo "Adding mismatches"
paf_add_mismatches -i ${working_dir}/combined.paf ${working_dir}/*.fa > ${working_dir}/mismatches.paf

# Run paf_trim
echo "Trimming"
paf_trim -i ${working_dir}/mismatches.paf -r 0.95 -t 0.0 > ${working_dir}/trimmed.paf

# Removing mismatches
echo "Removing mismatches"
paf_add_mismatches -i ${working_dir}/trimmed.paf -a > ${working_dir}/trimmed_no_mismatches.paf

# Run paf_chain
echo "Chaining"
paf_chain -i ${working_dir}/trimmed_no_mismatches.paf > ${working_dir}/chained.paf

# Run paf_tile
echo "Tiling"
paf_tile -i ${working_dir}/chained.paf > ${working_dir}/tiled.paf

# Get primary alignments
echo "Selecting primary alignments"
grep 'tp:A:P' ${working_dir}/tiled.paf > ${working_dir}/primary.paf

# Add back the mismatches
echo "Adding mismatches to the primary alignments"
paf_add_mismatches -i ${working_dir}/primary.paf ${working_dir}/*.fa > ${working_dir}/primary_mismatches.paf
