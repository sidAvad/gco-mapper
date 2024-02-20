#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name vcf2eig 
#SBATCH --output /home/sna53/gco_logs/vcf2eig.%a.log
#SBATCH --array=1-22

cd /home2/sna53/gene-conversions

counter=1
for chr in {1..22}
do
	if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
	then
		break 1
	fi
	counter=$(($counter+1))
done

source $HOME/miniconda3/bin/activate
python2 src/gdc/vcf2eigenstrat_phased.py -v rsync_data/admixedfile_phased_${chr}_200.vcf -o rsync_data/admixedfile_phased_${chr}_200

