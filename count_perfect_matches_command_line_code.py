#!/bin/bash
##
#SBATCH -o $logsdir/perfect_match_sh.${project}_%A_%a.log

source activate python_3.10.x

#Parameters
task=$SLURM_ARRAY_TASK_ID
sampleid=$(awk -F ',' -v task=$task 'FNR == task {print $1}' $intable)
ref=$(awk -F ',' -v task=$task 'FNR == task {print $2}' $intable)
fq=$(awk -F ',' -v task=$task 'FNR == task {print $3}' $intable)
outfile=$(awk -F ',' -v task=$task 'FNR == task {print $4}' $intable)

wdir=
refdir=$wdir/

echo $sampleid
echo $ref
echo $fq
echo $outfile

python $codedir/count_perfect_matches.py $sampleid $ref $fq $outfile $wdir $refdir 
