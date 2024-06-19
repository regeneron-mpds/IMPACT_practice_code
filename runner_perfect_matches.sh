#Fill in project here
project=


#This stays the same
wdir=$;codedir=$;intable=$;logsdir=$

# How to call this script:
sbatch --export=codedir=$codedir,intable=$intable,logsdir=$logsdir, -o $logsdir/perfect_match.${project}_%A_%a.log --array=1-20 $codedir/count_perfect_matches.sh