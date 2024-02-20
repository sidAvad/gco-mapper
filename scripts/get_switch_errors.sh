#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name get_switch_errors
#SBATCH --output /home/sna53/gco_logs/get_switch_error.%a.log
#SBATCH --array=1-22

cd /home2/sna53/gene-conversions

counter=1
for chr in {1..22}
do
    for size in 200
    do
	    if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
	    then
	    	break 2
	    fi
	    counter=$(($counter+1))
    done
done

VCFTOOLS_DIR=/programs/vcftools-v0.1.16/bin/vcftools
gunzip rsync_data/admixedfile_phased_${chr}_${size}.vcf.gz 

${VCFTOOLS_DIR} \
    --diff rsync_data/admixedfile_${chr}.vcf \
    --vcf rsync_data/admixedfile_phased_${chr}_${size}.vcf \
    --diff-switch-error \
    --out results/switch_errors_${chr}_${size}.txt
