#!/bin/bash
#SBATCH --partition regular,long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 32G
#SBATCH --exclude=cbsubscb18
#SBATCH --job-name simulate_YRICEU.%a.log
#SBATCH --output logs/simulate_YRICEU.%a.log

/programs/bin/labutils/mount_server cbsubscb18 /storage

SOURCE_DIR=v2_data/d00_inputs/source_data
#chrom="chr$SLURM_ARRAY_TASK_ID"

source /home/sna53/miniconda3/bin/activate popsim

SECONDS=0
python -m src.simulate_YRICEU 
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

echo '
import tskit
ts = tskit.load("simulation-source-stdpopsim-20k.trees")
print(ts.provenance(0))
' | python


