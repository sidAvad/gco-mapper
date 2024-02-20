#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name create_random_gcotables
#SBATCH --output /home/sna53/gco_logs/create_random_gcotables.%a.log
#SBATCH --array=1-10

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
gcotable=${GCO_DIR}/gcos_table_chr${chr}_hap${hap}.txt
output=${GCO_DIR}/gcos_table_chr${chr}_hap${hap}.rand.txt

source /home/sna53/miniconda3/bin/activate sciComp
python src/randomize_gcotable.py $gcotable $output
