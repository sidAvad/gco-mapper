#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks 45 
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name match_LAI_segments
#SBATCH --output logs/match_LAI_segments.hold.all.%a.log
#SBATCH --array=2

chr=1

hap=$SLURM_ARRAY_TASK_ID

source /home/sna53/miniconda3/bin/activate

admixed_file=rsync_data/gcos_all.recomb.gco_chr${chr}.hap${hap}.phgeno
CEU_file=rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.sim.phgeno
YRI_file=rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.sim.phgeno
ranges_file=testout.txt

#Filter the CEU and YRI files
#sh src/filter_phgeno.awk rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.hold.phgeno \
#    > rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.hold.flt.phgeno
#sh src/filter_phgeno.awk rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.hold.phgeno \
#    > rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.hold.flt.phgeno

YRI_file=rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.hold.phgeno
CEU_file=rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.hold.phgeno

window=1
line_number=0

echo $chr $hap $window

while read -r chrom hap start_phys end_phys start end gcos; do
	((line_number++))
	
    if [ "$line_number" -eq 1 ]; then
        # Skip the first line
        continue
    fi      
    
	echo "Raw input: $start|$end|$gcos"
    
    echo $start
    echo $end    
    
    sed -n "${start},${end}p" "$YRI_file" > tmp_data/YRI_holdout.$hap.$window.phgeno
    cat tmp_data/YRI_holdout.$hap.$window.phgeno | transpose > tmp_data/YRI_holdout.$hap.$window.trnsp
    sed -n "${start},${end}p" "$CEU_file" > tmp_data/CEU_holdout.$hap.$window.phgeno
    cat tmp_data/CEU_holdout.$hap.$window.phgeno | transpose > tmp_data/CEU_holdout.$hap.$window.trnsp
    
    cat tmp_data/YRI_holdout.$hap.$window.trnsp tmp_data/CEU_holdout.$hap.$window.trnsp \
        > tmp_data/ALL_holdout.$hap.$window.trnsp
    
    sed -n "${start},${end}p" "$admixed_file" > tmp_data/admixedfile.$hap.$window.phgeno
    cat tmp_data/admixedfile.$hap.$window.phgeno | transpose > tmp_data/admixedfile.$hap.$window.trnsp

    file1=tmp_data/ALL_holdout.$hap.$window.trnsp
    file2=tmp_data/admixedfile.$hap.$window.trnsp
    
	python src/match_LAIsegments.py -p $file1 -a $file2 -o \
        bestmatch_chr${chr}_hap${hap}_window${window}.phgeno > rsync_data/ndiffs_chr${chr}_hap${hap}_window${window}.hold.rand.txt &
    
    window=$(($window+1))

done < "$ranges_file"

# Wait for all background tasks to finish
wait
