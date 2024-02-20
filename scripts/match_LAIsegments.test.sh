#!/bin/bash

chr=1
hap=1
window_size=8000
GCO_DIR=v2_data/d00_inputs/gcos/${chr}
admixed_file=rsync_data/gcos_all.recomb_chr${chr}.hap${hap}.phgeno
CEU_sim_file=rsync_data/source_data/OutOfAfrica_4J17_CEU_chr${chr}.sim.phgeno
YRI_sim_file=rsync_data/source_data/OutOfAfrica_4J17_YRI_chr${chr}.sim.phgeno
ranges_file=rsync_data/gcos_all.recomb.gco_chr${chr}.rfmix.localanc.msp.tab_column_${hap}.split.txt

source /home/sna53/miniconda3/bin/activate

#file1=rsync_data/panel_test.txt
#file2=rsync_data/admixed_test.txt

head -n ${window_size} $CEU_sim_file > testCEU_${window_size}.phgeno
cat testCEU_${window_size}.phgeno | transpose > testCEU_${window_size}.trnsp
head -n ${window_size} $YRI_sim_file > testYRI_${window_size}.phgeno
cat testYRI_${window_size}.phgeno | transpose > testYRI_${window_size}.trnsp

cat testCEU_${window_size}.trnsp testYRI_${window_size}.trnsp > testALL_${window_size}.trnsp

head -n ${window_size} $admixed_file > testADMX_all_${window_size}.phgeno
cat testADMX_all_${window_size}.phgeno | transpose > testADMX_all_${window_size}.trnsp
head -n 1 testADMX_all_${window_size}.trnsp > testADMX_${window_size}.trnsp

python src/match_LAIsegments.py -p testALL_${window_size}.trnsp -a testADMX_${window_size}.trnsp -o output_${window_size}.txt > diffs_${window_size}.txt
