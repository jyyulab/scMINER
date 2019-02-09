#BSUB -P pbmc
#BSUB -M 24000
#BSUB -J SJARACNe
#BSUB -oo ZeroFiltered.out -eo ZeroFiltered.err
#BSUB -q pcgp

module load python/3.6.1

export PYTHON_PATH=/hpcf/apps/python/install/3.6.1/bin/python3.6
export scMINER_PATH=/research/projects/yu3grp/scRNASeq/yu3grp/TracyQian/scMINER/scMINER-master/

scMINER=/research/projects/yu3grp/scRNASeq/yu3grp/TracyQian/scMINER/scMINER-master/scMINER.py
indir=~/PBMCdemo/Sjaracne/

for i in $(ls -d */ | cut -f1 -d'/');do
python3 $scMINER SJARACNE ${i} $indir/${i}/*.exp $indir/${i}/tf/*.txt $indir/${i}/tf/ --c_threshold 0.01;
echo ${i} has been submitted!;
done


