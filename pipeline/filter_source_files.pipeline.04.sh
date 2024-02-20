#!/bin/bash
#SBATCH --partition regular,long7,long30
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name filer_source_files
#SBATCH --output logs/filter_source_files.log

chr=1

#Get alleles to be filtered out
#plink2 --vcf rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}.vcf \
#        --freq --out rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}

awk '{if ($5 > 0.05) print $2}' \
    rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr${chr}.frq > reference/snps_above_0.05_maf.backgroundgc.txt

#filter source vcfs
for pop in {'YRI','CEU'};
do
   /programs/bcftools-1.15.1-r/bin/bcftools view -i'ID=@reference/snps_above_0.05_maf.backgroundgc.txt' rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.vcf \
        > rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.filt.vcf 
    
    #Create eigenstrat format snp and phgeno files. 
    python2 src/gdc/vcf2eigenstrat_phased.py \
        -v rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.filt.vcf \
        -o rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.filt

    cat rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.filt.phgeno | transpose \
        > rsync_data/source_data/backgroundgc/OOA4J17_${pop}_chr${chr}.filt.phgeno.trnsp
done
