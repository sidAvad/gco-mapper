#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 64G
#SBATCH --job-name get_LDpatterns
#SBATCH --output logs/get_LDpatterns_perfectphase_rand.%a.log
#SBATCH --array=5

counter=1
for chr in 1
do
	for hap in {1..10}
	do
		if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
		then
			break 2
		fi
		counter=$(($counter+1))
	done
done


GCO_DIR=rsync_data/gco_tables
SOURCE_DIR=rsync_data/source_data
LD_RESULTS_DIR=results_rsync/LD_results_perfectphase_sim/${chr}
gcotable=${GCO_DIR}/gcos_table_chr${chr}_hap${hap}.rand.txt
ancestrytable=rsync_data/gcos_all.recomb.gco_chr${chr}.rfmix.localanc.msp.tab_column_${hap}.split.txt
admixed_phgeno=rsync_data/gcos_all.recomb.gco_chr${chr}.hap${hap}.phgeno
admixed_matched_phgeno=results_rsync/best_matches_perfectphase/bestmatch_chr${chr}_hap${hap}.phgeno
panelfile0=${SOURCE_DIR}/OutOfAfrica_4J17_YRI_chr${chr}.sim.phgeno
panelfile1=${SOURCE_DIR}/OutOfAfrica_4J17_CEU_chr${chr}.sim.phgeno
output_file=${LD_RESULTS_DIR}/gcos_table_chr${chr}_hap${hap}.rand.txt

#head -n 4950000 $admixed_phgeno > \
#    rsync_data/gcos_all.recomb.gco_chr${chr}.hap${hap}.phgeno.cropped
#head -n 4950000 $admixed_matched_phgeno > \
#   results_rsync/best_matches_imperfectphase/bestmatch_chr${chr}_hap${hap}.proc.phgeno.cropped
#head -n 4950000 $panelfile0 > \
#    ${SOURCE_DIR}/OutOfAfrica_4J17_YRI_chr${chr}.sim.phgeno.cropped
#head -n 4950000 $panelfile1 > \
#    ${SOURCE_DIR}/OutOfAfrica_4J17_CEU_chr${chr}.sim.phgeno.cropped

file1_lines=$(wc -l < "${panelfile0}.cropped")
file2_lines=$(wc -l < "${panelfile1}.cropped")
file3_lines=$(wc -l < "${admixed_matched_phgeno}.cropped")
file4_lines=$(wc -l < "${admixed_phgeno}.cropped")

if [ "$file1_lines" -eq "$file2_lines" ] && [ "$file1_lines" -eq "$file3_lines" ] && [ "$file1_lines" -eq "$file4_lines" ]; then
    echo "All files have the same number of lines: $file1_lines"
else
    echo "Number of lines in the files are not equal."
    echo "File 1: $file1_lines lines"
    echo "File 2: $file2_lines lines"
    echo "File 3: $file3_lines lines"
    echo "File 4: $file4_lines lines"
fi

source /home/sna53/miniconda3/bin/activate sciComp

#python v2_src/get_LDpatterns.py $gcotable $panelfile0 $panelfile1 ${admixed_matched_phgeno} ${admixed_phgeno} $output_file 
#python src/get_LDpatterns.py $gcotable $ancestrytable ${panelfile0}.cropped ${panelfile1}.cropped ${admixed_matched_phgeno}.cropped ${admixed_phgeno}.cropped ${output_file}
python src/get_LDpatterns_rand.py $gcotable $ancestrytable ${panelfile0}.cropped ${panelfile1}.cropped ${admixed_matched_phgeno}.cropped ${admixed_phgeno}.cropped ${output_file}
