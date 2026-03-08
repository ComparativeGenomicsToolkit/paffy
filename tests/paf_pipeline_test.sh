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
#wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simHuman.chr6 -O ${working_dir}/simHuman.chr6.fa
#wget https://raw.githubusercontent.com/UCSantaCruzComputationalGenomicsLab/cactusTestData/master/evolver/mammals/loci1/simMouse.chr6 -O ${working_dir}/simMouse.chr6.fa

# Run lastz
lastz ${working_dir}/mPanPan1_XY_1_5000000.fa[multiple][nameparse=darkspace] ${working_dir}/mPanTro3_XY_1_5000000.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz.paf
#lastz ${working_dir}/simHuman.chr6.fa[multiple][nameparse=darkspace] ${working_dir}/simMouse.chr6.fa[multiple][nameparse=darkspace] --step=4 --ambiguous=iupac,100,100 --ydrop=3000 --format=paf:minimap2 > ${working_dir}/lastz.paf

# Inverting so for every query:target alignment we have the mirror target:query alignment, this is used to
# make alignments symmetric for chaining
echo "Inverting"
paffy invert -i ${working_dir}/lastz.paf > ${working_dir}/inverted.paf

# Catting
echo "Catting forward and inverted pafs"
cat ${working_dir}/lastz.paf ${working_dir}/inverted.paf > ${working_dir}/combined.paf

# Split by query contig for parallel processing
echo "Splitting by query contig"
split_dir=${working_dir}/split
mkdir -p ${split_dir}
paffy split_file -q -i ${working_dir}/combined.paf -p ${split_dir}/

# List the parallel results
echo "Split results"
ls -l ${split_dir}

# Run add_mismatches -> chain -> tile -> trim in parallel per query contig
echo "Running add_mismatches -> chain -> tile -> trim in parallel"
pids=()
for split_file in ${split_dir}/*.paf; do
    base=$(basename ${split_file} .paf)
    out_dir=${working_dir}/parallel_${base}
    mkdir -p ${out_dir}
    (
        set -eo pipefail
        paffy add_mismatches -i ${split_file} ${working_dir}/*.fa \
          | paffy chain \
          | paffy tile \
          | paffy trim > ${out_dir}/trimmed.paf
    ) &
    pids+=($!)
done
# Wait for all parallel jobs and fail if any failed
for pid in "${pids[@]}"; do
    wait ${pid}
done

# Concatenate parallel results
echo "Concatenating parallel results"
cat ${working_dir}/parallel_*/trimmed.paf > ${working_dir}/trimmed.paf

# Report stats on the alignments after add_mismatches -> chain -> tile -> trim
echo "Reporting stats on trimmed alignments are as expected"
paffy view -i ${working_dir}/trimmed.paf ${working_dir}/*.fa -s -t

# Get primary alignments
echo "Selecting primary alignments"
paffy filter -i ${working_dir}/trimmed.paf -w 1  > ${working_dir}/primary.paf

# Report stats on the primary alignments before pruning
echo "Reporting stats on primary alignments are as expected"
paffy view -i ${working_dir}/primary.paf ${working_dir}/*.fa -s -t

# Run paffy chain again, to rechain given just the primary alignments
echo "Chaining primary alignments"
paffy chain -i ${working_dir}/primary.paf > ${working_dir}/primary_chained.paf

# Now filter out alignments in crappy primary chains
echo "Selecting primary alignments in good chains"
paffy filter -i ${working_dir}/primary_chained.paf -s 20000  > ${working_dir}/primary_final.paf

# Report stats on the primary alignments picked
echo "Reporting stats on primary alignments in good chains and check aligned bases and identity are as expected"
paffy view -i ${working_dir}/primary_final.paf ${working_dir}/*.fa -s -t -u 0.98 -v 16700000