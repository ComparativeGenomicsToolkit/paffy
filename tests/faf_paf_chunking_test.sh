#!/usr/bin/env bash

# Script to test running lastz using fasta chunking and paffy/faffy utils

# exit when any command fails
set -e

# log the commands its running
set -x

# Make a working directory
#working_dir=$(mktemp -d -t temp_chunking-XXXXXXXXXX)
working_dir=./temp_chunking
mkdir ${working_dir}

# Make sure we cleanup the temp dir
trap "rm -rf ${working_dir}" EXIT

# Get the sequences
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanPan1_XY_1_5000000.fa -O ${working_dir}/mPanPan1_XY_1_5000000.fa
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactusTestData/master/T2T_primate_PAR/mPanTro3_XY_1_5000000.fa -O ${working_dir}/mPanTro3_XY_1_5000000.fa

# Run faffy chunk on the two input sequences
echo "Chunking"
faffy chunk ${working_dir}/mPanPan1_XY_1_5000000.fa -c 1000000 -o 10000 -d ${working_dir}/temp_fastas_1
faffy chunk ${working_dir}/mPanTro3_XY_1_5000000.fa -c 1000000 -o 10000 -d ${working_dir}/temp_fastas_2

# Run lastz on all pairs of chunks
echo "Aligning"
for i in ${working_dir}/temp_fastas_1/*.fa
do
  for j in ${working_dir}/temp_fastas_2/*.fa
  do
    lastz ${i}[multiple][nameparse=darkspace] ${j}[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 >> ${working_dir}/lastz.paf
  done
done

#Now run lastz without chunking
lastz ${working_dir}/mPanPan1_XY_1_5000000.fa[multiple][nameparse=darkspace] ${working_dir}/mPanTro3_XY_1_5000000.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz_no_chunks.paf

# Use paffy dechunk to undo the effect of the chunking
echo "Dechunking"
paffy dechunk -i ${working_dir}/lastz.paf > ${working_dir}/lastz_dechunked.paf

# Remove dupes caused by overlaps
echo "Deduping"
paffy dedupe -i ${working_dir}/lastz_dechunked.paf > ${working_dir}/lastz_dechunked_dedupe.paf

# Report stats on the alignments
echo "Reporting stats on dechunked, deduped alignments and check aligned bases and identity are as expected"
paffy view -i ${working_dir}/lastz_dechunked_dedupe.paf ${working_dir}/*.fa -s -t -u 0.94 -v 22000000

# Report stats on the no chunk alignments
echo "Reporting stats on no chunk alignments for comparison"
paffy view -i ${working_dir}/lastz_no_chunks.paf ${working_dir}/*.fa -s -t -u 0.94 -v 22000000
