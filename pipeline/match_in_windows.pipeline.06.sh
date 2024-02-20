#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 128G
#SBATCH --job-name match_in_windows_hold
#SBATCH --output logs/match_in_windows.hold.%a.log
#SBATCH --array=1-5

counter=1
for window_size in {25000,50000,100000,200000,500000}
do
    if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
    then
        break
    fi
    counter=$((counter+1))
done

echo "window size =$window_size"

chr=1

YRIfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.hold.filt.phgeno"
CEUfile="rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr1.hold.filt.phgeno"
snpfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp"
admixedfile="rsync_data/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco.phgeno"
source /home/sna53/miniconda3/bin/activate py37


# Count lines in YRIfile
YRI_lines=$(wc -l < "$YRIfile")

# Count lines in CEUfile
CEU_lines=$(wc -l < "$CEUfile")

# Count lines in admixedfile
admixed_lines=$(wc -l < "$admixedfile")

# Check if all line counts are equal
if [ "$YRI_lines" -eq "$admixed_lines" ] && [ "$CEU_lines" -eq "$admixed_lines" ]; then
    echo "All files have equal line counts."
else
    echo "Line counts are not equal."
fi

SECONDS=0
python -m src.match_in_windows2 -a $admixedfile \
    -w $window_size \
    -Y $YRIfile \
    -C $CEUfile \
    -s $snpfile -o rsync_data/matches/300_individuals/matches_wbackgroundgc_windowsize${window_size}

duration=$SECOND
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
