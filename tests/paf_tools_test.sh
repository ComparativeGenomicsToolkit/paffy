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

# Run paffy left_align (moving gaps leaves identity and coverage as they were)
echo "paffy left_align minimum local alignment identity"
paffy left_align -i ${working_dir}/output.paf ${working_dir}/*.fa | paffy view ${working_dir}/*.fa -s -t -u 0.74 -v 530000

# Run paffy left_align on records carrying a tag paffy does not parse: everything but the cigar comes through as it was
echo "paffy left_align keeps the rest of each record"
awk 'BEGIN{OFS="\t"} {print $0, "rc:Z:x"NR}' ${working_dir}/output.paf > ${working_dir}/tagged.paf
paffy left_align -i ${working_dir}/tagged.paf ${working_dir}/*.fa > ${working_dir}/tagged.left.paf
diff <(sed 's/\tcg:Z:[^\t]*//' ${working_dir}/tagged.paf) <(sed 's/\tcg:Z:[^\t]*//' ${working_dir}/tagged.left.paf)

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

# Run paffy unanchor on a toy chromosome: a (CA)n locus on the reference path CHM13#0#chr1 = s1 s2 s3, which two
# haplotypes (one on each strand) read through a longer alt allele s4. Two compensating contigs pass the gate, so
# the locus's hub intervals are cleared from every record, the reference's included (expected outputs are those
# of the anchor-trim pilot's Python reference)
echo "paffy unanchor toy"
tr ' ' '\t' > ${working_dir}/toy.gfa <<'TOY'
S s1 * LN:i:100 SN:Z:CHM13#0#chr1 SO:i:0 SR:i:0
S s2 * LN:i:30 SN:Z:CHM13#0#chr1 SO:i:100 SR:i:0
S s3 * LN:i:100 SN:Z:CHM13#0#chr1 SO:i:130 SR:i:0
S s4 * LN:i:40 SN:Z:H1#1#ctg1 SO:i:100 SR:i:1
TOY
cat > ${working_dir}/toy.fa <<'TOY'
>s1
GCTAAAGACAATTACATAACATACACGTCAGCACGAAACTTGTTGGCCCAGTGTGAATCGCTTAAGGGTTAAGTAAGTGTGATGCATACGCCTTTACTTG
>s2
CACACACACACACACACACACACACACACA
>s3
CTGTGTCCACCCCATCGGACTGGCATTTTTATTACACTCAGAAACAGAACTCGGGTAATTTTGACAGGTCACGCAGAGGCGCGCCCTCCTGAAGTGCGTG
>s4
CACACACACACACACACACACACACACACACACACACACA
TOY
tr ' ' '\t' > ${working_dir}/toy.paf <<'TOY'
id=CHM13|chr1 230 0 100 + id=_MINIGRAPH_|s1 100 0 100 100 100 60 tp:A:P cg:Z:100=
id=CHM13|chr1 230 100 130 + id=_MINIGRAPH_|s2 30 0 30 30 30 60 tp:A:P cg:Z:30=
id=CHM13|chr1 230 130 230 + id=_MINIGRAPH_|s3 100 0 100 100 100 60 tp:A:P cg:Z:100=
id=H1.1|ctg1 240 0 100 + id=_MINIGRAPH_|s1 100 0 100 75 100 60 NM:i:0 AS:i:100 tp:A:P cg:Z:100M rc:Z:x
id=H1.1|ctg1 240 100 140 + id=_MINIGRAPH_|s4 40 0 40 40 40 60 tp:A:P cg:Z:40=
id=H1.1|ctg1 240 140 240 + id=_MINIGRAPH_|s3 100 0 99 98 101 60 tp:A:P cg:Z:5=1D45=2I48= cs:Z:junk
id=H1.1|ctg1 240 100 130 + id=_MINIGRAPH_|s2 30 0 30 30 30 0 tp:A:S cg:Z:30=
id=H2.1|ctg1 240 0 100 - id=_MINIGRAPH_|s3 100 0 100 100 100 60 tp:A:P cg:Z:100=
id=H2.1|ctg1 240 100 140 - id=_MINIGRAPH_|s4 40 0 40 40 40 60 tp:A:P cg:Z:40=
id=H2.1|ctg1 240 140 240 - id=_MINIGRAPH_|s1 100 0 100 100 100 60 tp:A:P cg:Z:60=1X39=
TOY
toy_args="-i ${working_dir}/toy.paf -n ${working_dir}/toy.fa -g ${working_dir}/toy.gfa -r CHM13"
paffy unanchor ${toy_args} -o ${working_dir}/toy.cut.paf -L ${working_dir}/toy.loci.bed -b ${working_dir}/toy.hub.bed \
    -v ${working_dir}/toy.veto.bed -q ${working_dir}/toy.query.bed
diff ${working_dir}/toy.loci.bed - <<'TOY'
#contig	start	end	locus	period	class	boundary	nHaps	compHaps	insideDeps	compAltnode	compReplace	compCigar	minString	maxString	refIntervals	altNodes	alleles	contigEnds
chr1	81	150	0	2	bubble	1	2	2	3	2	0	0	75	85	3	1	78x1,79x1	0
TOY
diff ${working_dir}/toy.hub.bed - <<'TOY'
#target	start	end	locus	kind	class	boundary	nHaps	compHaps	maxString
id=_MINIGRAPH_|s1	81	100	0	ref	bubble	1	2	2	85
id=_MINIGRAPH_|s2	0	30	0	ref	bubble	1	2	2	85
id=_MINIGRAPH_|s3	0	20	0	ref	bubble	1	2	2	85
id=_MINIGRAPH_|s4	0	40	0	alt	bubble	1	2	2	85
TOY
# the cut: fragments on both strands keep every tag but NM AS cs (and the rest of the record) byte for byte; M
# columns are rescaled (81 * 75/100 = 60.75 -> 61); gaps next to a cut go; the secondary on s2 is dropped
diff ${working_dir}/toy.cut.paf - <<'TOY'
id=CHM13|chr1	230	0	81	+	id=_MINIGRAPH_|s1	100	0	81	81	81	60	tp:A:P	cg:Z:81=
id=CHM13|chr1	230	150	230	+	id=_MINIGRAPH_|s3	100	20	100	80	80	60	tp:A:P	cg:Z:80=
id=H1.1|ctg1	240	0	81	+	id=_MINIGRAPH_|s1	100	0	81	61	81	60	tp:A:P	cg:Z:81M	rc:Z:x
id=H1.1|ctg1	240	159	240	+	id=_MINIGRAPH_|s3	100	20	99	79	81	60	tp:A:P	cg:Z:31=2I48=
id=H2.1|ctg1	240	0	80	-	id=_MINIGRAPH_|s3	100	20	100	80	80	60	tp:A:P	cg:Z:80=
id=H2.1|ctg1	240	159	240	-	id=_MINIGRAPH_|s1	100	0	81	80	81	60	tp:A:P	cg:Z:60=1X20=
TOY
# the query bases that lost every hub anchor, in each contig's own (forward) coordinates
diff ${working_dir}/toy.query.bed - <<'TOY'
id=CHM13|chr1	81	150
id=H1.1|ctg1	81	159
id=H2.1|ctg1	80	159
TOY
[ "$(grep -vc '^#' ${working_dir}/toy.veto.bed)" = 0 ]
# the same with 3 threads, and a second run over the cut changes nothing (the reference's own anchors are gone)
paffy unanchor ${toy_args} -t 3 | diff ${working_dir}/toy.cut.paf -
paffy unanchor -i ${working_dir}/toy.cut.paf -n ${working_dir}/toy.fa -g ${working_dir}/toy.gfa -r CHM13 | diff ${working_dir}/toy.cut.paf -
# --dryRun copies the PAF
paffy unanchor ${toy_args} -d -q ${working_dir}/toy.dry.bed | diff ${working_dir}/toy.paf -
[ ! -s ${working_dir}/toy.dry.bed ]
# no reference contig (a chrOther-like bin): the PAF is copied unchanged and the run succeeds
paffy unanchor ${toy_args/CHM13/GRCh38} | diff ${working_dir}/toy.paf -
# inconsistent inputs exit 2: a target that is not the hub, and a node whose FASTA and GFA lengths disagree
sed 's/id=_MINIGRAPH_|s4/id=X|s4/' ${working_dir}/toy.paf > ${working_dir}/toy.bad.paf
set +e
paffy unanchor -i ${working_dir}/toy.bad.paf -n ${working_dir}/toy.fa -g ${working_dir}/toy.gfa -r CHM13 > /dev/null 2>&1
bad_target=$?
sed 's/LN:i:40/LN:i:41/' ${working_dir}/toy.gfa > ${working_dir}/toy.bad.gfa
paffy unanchor -i ${working_dir}/toy.paf -n ${working_dir}/toy.fa -g ${working_dir}/toy.bad.gfa -r CHM13 > /dev/null 2>&1
bad_length=$?
set -e
[ ${bad_target} = 2 ] && [ ${bad_length} = 2 ]
