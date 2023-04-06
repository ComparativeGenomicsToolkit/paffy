#!/usr/bin/env bash

# Script to test pafs

# exit when any command fails
set -e

# log the commands its running
set -x

# Make a working directory
#working_dir=$(mktemp -d -t temp_chains-XXXXXXXXXX)
working_dir=./temp_chains
mkdir -p ${working_dir}

# Make sure we cleanup the temp dir
trap "rm -rf ${working_dir}" EXIT

# Get the sequences
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanPan1_XY_1_5000000.fa -O ${working_dir}/mPanPan1_XY_1_5000000.fa
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanTro3_XY_1_5000000.fa -O ${working_dir}/mPanTro3_XY_1_5000000.fa
#wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simCow.chr6 -O ${working_dir}/simCow.chr6.fa
#wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simDog.chr6 -O ${working_dir}/simDog.chr6.fa

# Run lastz
lastz ${working_dir}/mPanPan1_XY_1_5000000.fa[multiple][nameparse=darkspace] ${working_dir}/mPanTro3_XY_1_5000000.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz.paf
#lastz ${working_dir}/simCow.chr6.fa[multiple][nameparse=darkspace] ${working_dir}/simDog.chr6.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz.paf

# Inverting so for every query:target alignment we have the mirror target:query alignment, this is used to
# make alignments symmetric for chaining
echo "Inverting"
paffy invert -i ${working_dir}/lastz.paf > ${working_dir}/inverted.paf

# Catting
echo "Catting forward and inverted pafs"
cat ${working_dir}/lastz.paf ${working_dir}/inverted.paf > ${working_dir}/combined.paf

# Run paffy add_mismatches
echo "Adding mismatches"
paffy add_mismatches -i ${working_dir}/combined.paf ${working_dir}/*.fa > ${working_dir}/mismatches.paf

# Run paffy chain
echo "Chaining"
paffy chain -i ${working_dir}/mismatches.paf > ${working_dir}/chained.paf

# Run paffy tile
echo "Tiling"
paffy tile -i ${working_dir}/chained.paf > ${working_dir}/tiled.paf

# Run paffy trim
echo "Trimming"
paffy trim -i ${working_dir}/tiled.paf > ${working_dir}/trimmed.paf

# Get primary alignments
echo "Selecting primary alignments"
grep 'tp:A:P' ${working_dir}/trimmed.paf > ${working_dir}/primary.paf

# Report stats on the primary alignments picked
echo "Reporting stats on primary alignments and check aligned bases and identity are as expected"
paffy view -i ${working_dir}/primary.paf ${working_dir}/*.fa -s -t -u 0.98 -v 17400000