#!/bin/bash
#SBATCH --partition regular,long7,long30
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 32G
#SBATCH --job-name sample_admixed_genotypes
#SBATCH --output logs/sample_admixed_genotypes.log

chr=1

YRIfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.phgeno.trnsp"
CEUfile="rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr1.sim.filt.phgeno.trnsp"
snpfile="rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp"
outprefix="rsync_data/simulated_data/100_individuals_testgc/admixed_w_backgroundgc"
gcofilelist="rsync_data/bpfiles/gcofilelist.txt"
bpfile_all="rsync_data/bpfiles/gcos_all.recomb_chr1.mrk.bp"

#Remove composite files and results from previous runs
rm $gcofilelist
rm $bpfile_all
rm $outprefix*

#Create composite files
for ind in {1..100}
do
    cat rsync_data/bpfiles/gcos_${ind}_0.recomb_chr1.mrk.bp
    cat rsync_data/bpfiles/gcos_${ind}_1.recomb_chr1.mrk.bp
    echo ""
done > $bpfile_all

for ind in {1..100} #Remove gcos with no markers in them
do
    grep -v "\[\]" rsync_data/bpfiles/gcos_${ind}_0.gco_chr1.mrk.bp > rsync_data/bpfiles/gcos_${ind}_0.gco_chr1.mrk.onlysnps.bp 
    grep -v "\[\]" rsync_data/bpfiles/gcos_${ind}_1.gco_chr1.mrk.bp > rsync_data/bpfiles/gcos_${ind}_1.gco_chr1.mrk.onlysnps.bp 
done

for ind in {1..100}
do
    echo "rsync_data/bpfiles/gcos_${ind}_0.gco_chr1.mrk.onlysnps.bp"
    echo "rsync_data/bpfiles/gcos_${ind}_1.gco_chr1.mrk.onlysnps.bp"
done > $gcofilelist

source /home/sna53/miniconda3/bin/activate sciComp
python -m src.sample_admixed_genotypes -b $bpfile_all -g $gcofilelist -y $YRIfile -c $CEUfile -s $snpfile -o $outprefix

cat ${outprefix}_recomb_gco.txt | transpose > ${outprefix}_recomb_gco.phgeno
