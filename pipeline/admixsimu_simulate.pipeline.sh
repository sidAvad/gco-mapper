#!/bin/bash
#SBATCH --partition regular,long7,long30
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name admixSimu_simulate
#SBATCH --output logs/admixSimu_simulate.%a.log
#SBATCH --array=1

#/programs/bin/labutils/mount_server cbsubscb18 /storage
#cd /fs/cbsubscb18/storage/siddharth/gene-conversions


chr=$SLURM_ARRAY_TASK_ID

#put this in if the first line of the .asu.mrk.bp file is not 'YRI CEU\n'
#sed -i '1s/^/YRI CEU\n/' rsync_data/bpfiles/gcos_all.recomb_chr${chr}.asu.mrk.bp

src/simu-mix.pl -bp rsync_data/bpfiles/gcos_all.recomb_chr${chr}.asu.mrk.bp \
    rsync_data/simulated_data/simulated.recomb.chr${chr} \
    -CEU rsync_data/source_data/OutOfAfrica_4J17_CEU_chr1.sim.filt.phgeno \
	-YRI rsync_data/source_data/OutOfAfrica_4J17_YRI_chr1.sim.filt.phgeno

