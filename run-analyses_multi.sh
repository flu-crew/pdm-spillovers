# Commands for pdm09-spillover analyses. The script is designed to be run in parallel multiple times to obtain multiple replicates.

# Construct the tree for human+swine seqs:
analysis="PDM-HA-USA-09-21-Mar2022-merged-cds"
aln="trees/pdm09_US_all_test.aln"
run=$1  # Replicate number (1-20).

cat human-swine-analysis/data/swine-$analysis.aln human-swine-analysis/data/human-$analysis.aln > $aln
directory="trees/results/$analysis-$run/"
mkdir $directory
aln2=${directory}pdm09_US_all.aln
cp $aln $aln2
aln=$aln2

./run-iqtree.py $aln -m GTR+F+R5 -pre ${directory}iqtree
iqtree_tree=${directory}iqtree.treefile
#fasttree -nt -gtr -gamma $aln > $iqtree_tree  #  A faster alternative to IQ-Tree if time/resources are limited.


## Root the tree:
./treetime-root.py $iqtree_tree $aln
rooted_tree=$aln.rooted.tre
dates=$aln.dates.csv

## Infer time-scaled tree:
outdir="${directory}pdm09_US_all_timetree"
treetime --tree $rooted_tree --aln $aln --dates $dates --outdir $outdir
time_tree=$aln.timetree.tre
python dendropy-format-converter.py $outdir/timetree.nexus nexus newick
cp $outdir/timetree.newick $time_tree

## Infer ancestral hosts (mugration model):
outdir="${directory}pdm09_US_all_mugration"
treetime mugration --tree $time_tree --states trees/$analysis.hosts.csv --attribute host --outdir $outdir
host_tree="${directory}pdm09_US_all.timetree.hosts.tre"
cp $outdir/annotated_tree.nexus $host_tree

## Infer HA1 ancestral substitutuons
ha1_aln="${directory}pdm09_US_all.ha1.aln"
flutile trim ha1 --conversion dna2aa --subtype H1 $aln > $ha1_aln
outdir="${directory}pdm09_US_all_aasub_ancestral"
treetime ancestral --tree $time_tree --aln $ha1_aln --aa --outdir $outdir
