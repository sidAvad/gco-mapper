#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --mem-per-cpu 16G
#SBATCH --job-name filter_snps
#SBATCH --output logs/filter_snps.all.log

chr=1
typ=sim

#First the YRI file is filtered out and line numbers to be deleted are recorded
#Next the admixed file and CEU file is filtered out using the recorded line-numbers
#Importantly, the snp file is also filtered using the recorded line-numebers; \
    #it can be used to translate physical to genetic position

echo $chr $typ

CEU_file=rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.$typ.phgeno
YRI_file=rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.$typ.phgeno

#Get line numbers to be filtered out of YRI file
sh src/filter_phgeno_linenumbers.awk rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.$typ.phgeno \
    > rsync_data/source_data/filtered_lines.sed

echo "lines to be deleted are recorded"

#Filter out those lines from CEU and admixed file
sed -i -f rsync_data/source_data/filtered_lines.sed \
    $YRI_file > rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.$typ.flt.phgeno
sed -i -f rsync_data/source_data/filtered_lines.sed \
    $CEU_file > rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.$typ.flt.phgeno

#Filter out lines in snp file (important because this snp file is \
    #required to translate back and forth between physical and marker positions)
sed -i -f rsync_data/source_data/filtered_lines.sed \
    rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.snp > rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.flt.snp

# Get the line count of each file
line_count1=$(wc -l < rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.$typ.flt.phgeno)
line_count2=$(wc -l < rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.$typ.flt.phgeno)
line_count3=$(wc -l < rsync_data/gcos_all.recomb.gco_chr${chr}.hap${hap}.flt.phgeno)
line_count4=$(wc -l < rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.flt.snp)

# Compare the line counts
if [[ $line_count1 -eq $line_count2 && $line_count2 -eq $line_count3 && $line_count3 -eq $line_count4 ]]; then
    echo "The line numbers in the three files are equal."
else
    echo "The line numbers in the three files are not equal."
fi
