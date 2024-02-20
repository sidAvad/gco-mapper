#!/bin/bash
#SBATCH --partition regular
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --time 1-0:00:00
#SBATCH --exclude=cbsubscb18,cbsubscb14
#SBATCH --job-name gco-sim
#SBATCH --output logs/gco-sim-%a.log
#SBATCH --array=1-1000

/programs/bin/labutils/mount_server cbsubscb18 /storage

hap=$SLURM_ARRAY_TASK_ID
python gcosim/gcosim.py 7 ${hap} gcos reference/refined_mf.simmap
#
mv gcos_${hap}.* rsync_data/bpfiles/
#
python src/splitChrom_AdmixSimu.py -i rsync_data/bpfiles/gcos_${hap}.recomb.bp
python src/phys2mark_AdmixSimu.py -i rsync_data/bpfiles/gcos_${hap}.recomb_chr1.bp \
    -s rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp
head -n 1 rsync_data/bpfiles/gcos_${hap}.recomb_chr1.mrk.bp > rsync_data/bpfiles/gcos_${hap}_0.recomb_chr1.mrk.bp
tail -n 1 rsync_data/bpfiles/gcos_${hap}.recomb_chr1.mrk.bp > rsync_data/bpfiles/gcos_${hap}_1.recomb_chr1.mrk.bp

#gco file
awk '($5==0 || NR==1){print $0}' rsync_data/bpfiles/gcos_${hap}.gco.bp >  rsync_data/bpfiles/gcos_${hap}_0.gco.bp
awk '($5==1 || NR==1){print $0}' rsync_data/bpfiles/gcos_${hap}.gco.bp >  rsync_data/bpfiles/gcos_${hap}_1.gco.bp
awk '($1==1 || NR==1){print $0}' rsync_data/bpfiles/gcos_${hap}_0.gco.bp >  rsync_data/bpfiles/gcos_${hap}_0.gco_chr1.bp
awk '($1==1 || NR==1){print $0}' rsync_data/bpfiles/gcos_${hap}_1.gco.bp >  rsync_data/bpfiles/gcos_${hap}_1.gco_chr1.bp

python src/phys2mark_AdmixSimu.py -i rsync_data/bpfiles/gcos_${hap}_0.gco_chr1.bp \
    -s rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp -g
python src/phys2mark_AdmixSimu.py -i rsync_data/bpfiles/gcos_${hap}_1.gco_chr1.bp \
    -s rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp -g

#python src/phys2mark_AdmixSimu.py -i rsync_data/bpfiles/gcos_${hap}_0.gco_chr1.mrk.bp \
#    -s rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp -g
