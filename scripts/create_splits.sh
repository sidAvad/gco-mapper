#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name create_splits
#SBATCH --output logs/create_splits.%a.log
#SBATCH --array=1-2

/programs/bin/labutils/mount_server cbsubscb18 /storage
cd /fs/cbsubscb18/storage/siddharth/gene-conversions

counter=1
for chr in {1..22}
do
	for group in {"sim","hold"}
	do
		if [[ $counter -eq $SLURM_ARRAY_TASK_ID ]];then
			break 2
		fi
		counter=$(($counter+1))
	done
done


SOURCE_DIR=rsync_data/source_data/backgroundgc

#activate conda
source $HOME/miniconda3/bin/activate


echo "chr=${chr} & group=${group}"

#subset vcfs to simulation and hold sets (full sets are created in subset_YRICEU.sh)
/programs/bcftools-1.15.1-r/bin/bcftools view -S reference/${group}_YRI.ids.txt ${SOURCE_DIR}/OOA4J17_YRI_chr${chr}.filt.vcf \
    -o ${SOURCE_DIR}/OOA4J17_YRI_chr${chr}.${group}.filt.vcf
/programs/bcftools-1.15.1-r/bin/bcftools view -S reference/${group}_CEU.ids.txt ${SOURCE_DIR}/OOA4J17_CEU_chr${chr}.filt.vcf \
    -o ${SOURCE_DIR}/OOA4J17_CEU_chr${chr}.${group}.filt.vcf


python2 ~/siddharth/Tools/python-packages_iain/gdc/vcf2eigenstrat_phased.py \
	-v ${SOURCE_DIR}/OOA4J17_YRI_chr${chr}.${group}.filt.vcf \
	-o ${SOURCE_DIR}/OOA4J17_YRI_chr${chr}.${group}.filt

python2 ~/siddharth/Tools/python-packages_iain/gdc/vcf2eigenstrat_phased.py \
	-v ${SOURCE_DIR}/OOA4J17_CEU_chr${chr}.${group}.filt.vcf \
	-o ${SOURCE_DIR}/OOA4J17_CEU_chr${chr}.${group}.filt
