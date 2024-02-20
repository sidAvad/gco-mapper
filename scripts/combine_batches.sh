#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name collapse_matched_segments 
#SBATCH --output /home/sna53/gco_logs/collapse_matched_segments_imperfectphase.%a.log
#SBATCH --array=1-10

chr=1
counter=1
for hap in {1..10}
do
    if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
    then
    	break 2
    fi
    counter=$(($counter+1))
done

#Collapse phgenos across windows to make a single file
output_file="results_rsync/best_matches_imperfectphase/bestmatch_chr${chr}_hap${hap}.phgeno"
> $output_file

for batch in {1..15}
do
    for window in {1..60}
    do
        echo "batch : ${batch}"
        echo "window : ${window}"
        input_file="results_rsync/best_matches_imperfectphase/bestmatch_chr${chr}_hap${hap}_batch${batch}_window${window}.phgeno"
        
		if grep -qv '^.$' "$filename"; then
            echo "file $filename does not have exactly one character per line"
        fi
		
        # Check if the input file exists before concatenating
        if [ -f "$input_file" ]; then
            echo "Concatenating $input_file"
            cat "$input_file" >> "$output_file"
        fi	
		
    done
done
