#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 16G
#SBATCH --exclude=cbsubscb18
#SBATCH --job-name collate_results
#SBATCH --output logs/collate_results.%a.log
#SBATCH --array=1-5

counter=1
for window_size in {25000,50000,100000,200000,500000}
do
    if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
    then
        break
    fi
    counter=$((counter+1))
done
#window_size=100000
source /home/sna53/miniconda3/bin/activate sciComp
python src/collate_results.py 200 $window_size
