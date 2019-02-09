#BSUB -P demo
#BSUB -n 1
#BSUB -M 12000
#BSUB -oo demo.out -eo demo.err
#BSUB -J MICA
#BSUB -q pcgp

python3=/hpcf/apps/python/install/3.6.1/bin/python3.6
scMINER=/research/projects/yu3grp/scRNASeq/yu3grp/TracyQian/scMINER/scMINER-master/scMINER.py

sample="PBMC_Demo"
species="hg19"
indir=~/PBMC_Demo/

$python3 $scMINER MIE Pipeline $sample $indir/PBMC_Demo_MICA_input_mini.txt $indir $sample

$python3 $scMINER MICA Clust $sample $indir/scMINER_MIE_out/PBMC_Demo.whole.h5 $indir/scMINER_PBMC_Demo/scMINER_MIE_out/PBMC_Demo_mi.h5 $indir/scMINER_MICA/ $sample --k 3 4 5 6 --perplexity 30 --retransformation False
