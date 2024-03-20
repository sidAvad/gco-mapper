#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 64G
#SBATCH --exclude=cbsubscb18
#SBATCH --job-name sub_match
#SBATCH --output logs/sub_match.log

#YRIfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.hold.filt.phgeno"
#CEUfile="rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr1.hold.filt.phgeno"
#snpfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp"
#admixedfile="rsync_data/simulated_data/100_individuals_testgc/admixed_w_backgroundgc_recomb_gco.phgeno"
#
#python -m src.sub_match -W results/100_individuals_testgc/results_windowsize200000.txt \
#    -w 100000 \
#    -a $admixedfile \
#    -Y $YRIfile
#    -C $CEUfile \
#    -s $snpfile \
#    -o testout.txt

/home/sna53/miniconda3/bin/activate py37
python -m unittest src/test_sub_match.py

