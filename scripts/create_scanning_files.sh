#!/bin/bash

################

# Run the python script to turn gco_tables into scanning files by creating windows around 
# gcos ( with physical position and marker numbers ) as well as random windows. 
# This will allow quicker testing of matching algorithm by focusing on these regions. 

################

source /home/sna53/miniconda3/bin/activate sciComp

GCO_DIR=rsync_data/gco_tables
gcotable=${GCO_DIR}/gcos_table_chr${chr}_hap${hap}.txt
mapfile=reference/sex-averaged.simmap
output=testout.txt

python src/get_scanning_file.py -i $gco_table -m $mapfile -o $output


