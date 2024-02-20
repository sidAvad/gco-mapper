#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name process_files 
#SBATCH --output /home/sna53/gco_logs/process_files.%a.log
#SBATCH --array=1-10

# Function to process lines and split multi-character lines
process_file() {
    input_file="$1"
    output_file="$2"
    
    while IFS= read -r line; do
        if [ "${#line}" -gt 1 ]; then
            for ((i = 0; i < ${#line}; i++)); do
                echo "${line:$i:1}" >> "$output_file"
            done
        else
            echo "$line" >> "$output_file"
        fi
    done < "$input_file"
}

# Loop over values of hap from 1 to 10
hap=$SLURM_ARRAY_TASK_ID
echo ${hap}
input_file="results_rsync/best_matches_imperfectphase/bestmatch_chr1_hap${hap}.phgeno"
output_file="results_rsync/best_matches_imperfectphase/bestmatch_chr1_hap${hap}.proc.phgeno"
> $output_file
process_file "$input_file" "$output_file"

echo "Processing completed."

