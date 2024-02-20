#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 32G
#SBATCH --exclude=cbsubscb18
#SBATCH --job-name subset_YRICEU
#SBATCH --output logs/subset_YRICEU.%a.log
#SBATCH --array=1

chr=$SLURM_ARRAY_TASK_ID

#source /home/sna53/miniconda3/bin/activate popsim

#echo '
#import tskit
#ts = tskit.load("simulation-source-stdpopsim-20k.trees")
#pop_id_to_name = {pop.id : pop.metadata["name"] for pop in ts.populations()}
#print(pop_id_to_name)
#
#indiv_pop = [pop_id_to_name[ind.population] for ind in ts.individuals()]
#print(indiv_pop)
#
#indiv_name = [pop + str(i) for i, pop in enumerate(indiv_pop)]
#print(indiv_name)
#ts.write_vcf(open("OOA4J17_withgc_chr'${chr}'.vcf", "w"),individual_names=indiv_name)
#' | python

/programs/bcftools-1.15.1-r/bin/bcftools view -S reference/OOA4J17_YRI.samples.txt \
	rsync_data/source_data/backgroundgc/OOA4J17_withgc_chr${chr}.vcf \
	> rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}.vcf
/programs/bcftools-1.15.1-r/bin/bcftools view -S reference/OOA4J17_CEU.samples.txt \
	rsync_data/source_data/backgroundgc/OOA4J17_withgc_chr${chr}.vcf \
	> rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr${chr}.vcf

python2 src/gdc/vcf2eigenstrat_phased.py \
	-v rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}.vcf \
	-o rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}
python2 ~/siddharth/Tools/python-packages_iain/gdc/vcf2eigenstrat_phased.py \
	-v rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr${chr}.vcf \
	-o rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr${chr}
