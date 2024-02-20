#!/bin/bash
#SBATCH --partition regular,long7
#SBATCH --nodes 1
#SBATCH --ntasks 16
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name phase_gcofiles
#SBATCH --output logs/phase_gcofiles.backgroundgc.log

chr=1
BCFTOOLS_PATH=/programs/bcftools-1.15.1-r/bin/bcftools
TABIX_PATH=/programs/htslib-1.16/bin/tabix 
RSYNCDIR=/home2/sna53/gene-conversions/rsync_data
YRI_samples=reference/hold_YRI.ids.txt
CEU_samples=reference/hold_CEU.ids.txt
size=5000

source /home/sna53/miniconda3/bin/activate sciComp

#unphase admixed file (convert phgeno to geno)
python src/eig2eig.py $RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco.phgeno

#Create vcf from unphased admixed file
python src/eig2vcf.py -r $RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco \
    -s $RSYNCDIR/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt \
    > $RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco.vcf

half_size=$((size / 2))

## Randomly select samples from CEU and YRI groups
selected_yri_samples=$(cat "$YRI_samples" | shuf -n "$half_size" | tr '\n' ',')
selected_ceu_samples=$(cat "$CEU_samples" | shuf -n "$half_size" | tr '\n' ',')

#Select the samples and write to VCF file
"$BCFTOOLS_PATH" view -s "$selected_yri_samples" \
    --force-samples "${RSYNCDIR}/source_data/backgroundgc/OOA4J17_YRI_chr1.hold.filt.vcf" \
    -O v -o "${RSYNCDIR}/source_data/backgroundgc/YRI_phasingreference_${chr}_${size}.vcf"
"$BCFTOOLS_PATH" view -s "$selected_ceu_samples" \
    --force-samples "${RSYNCDIR}/source_data/backgroundgc/OOA4J17_CEU_chr1.hold.filt.vcf" \
    -O v -o "${RSYNCDIR}/source_data/backgroundgc/CEU_phasingreference_${chr}_${size}.vcf"

#Combined vcfs to create reference file
"$BCFTOOLS_PATH" merge ${RSYNCDIR}/source_data/backgroundgc/YRI_phasingreference_${chr}_${size}.filt.vcf.gz \
    ${RSYNCDIR}/source_data/backgroundgc/CEU_phasingreference_${chr}_${size}.filt.vcf.gz \
    > ${RSYNCDIR}/source_data/backgroundgc/YRICEU_phasingreference_${chr}_${size}.filt.vcf

file1=$RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco.vcf
file2=$RSYNCDIR/source_data/backgroundgc/YRICEU_phasingreference_${chr}_${size}.filt.vcf
wc1=$(grep -v '^#' $file1 | wc -l)
wc2=$(grep -v '^#' $file2 | wc -l)
echo "nlines file 1 = $wc1 and nlines file 2 = $wc2"
[[ $wc1 -eq $wc2 ]] \
    || { echo "The number of lines is different"; exit 1; }

java -jar beagle/beagle.22Jul22.46e.jar \
    gt=$file1 \
    ref=$file2 \
    out=$RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco_phased_${chr}_${size} \
    map=$RSYNCDIR/reference/plink.GRCh37.map/plink.chr${chr}.GRCh37.map chrom=${chr} nthreads=16

source $HOME/miniconda3/bin/activate
python2 src/gdc/vcf2eigenstrat_phased.py \
    -v $RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco_phased_${chr}_${size}.vcf.gz \
    -o $RSYNCDIR/simulated_data/300_individuals/admixed_w_backgroundgc_recomb_gco_phased_${chr}_${size}.phgeno

#Some code for bgzip and tabix
#for file in {"${RSYNCDIR}/source_data/no_backgroundgc/YRI_phasingreference_${chr}_${size}.vcf","${RSYNCDIR}/source_data/no_backgroundgc/YRI_phasingreference_${chr}_${size}.vcf"}
#do
#    echo $file
#    /programs/htslib-1.11/bin/bgzip -c $file \
#        > $file.gz
#    /programs/htslib-1.11/bin/tabix -f -p vcf \
#        $file.gz
#done
#
#file=${RSYNCDIR}/source_data/backgroundgc/YRI_phasingreference_${chr}_${size}
#/programs/bcftools-1.15.1-r/bin/bcftools view --max-alleles 2 --exclude-types indels $file.vcf > $file.filt.vcf
#
#echo $file
#/programs/htslib-1.11/bin/bgzip -c $file.filt.vcf \
#    > $file.filt.vcf.gz
#/programs/htslib-1.11/bin/tabix -f -p vcf \
#    $file.filt.vcf.gz
#
#file=${RSYNCDIR}/source_data/backgroundgc/CEU_phasingreference_${chr}_${size}
#/programs/bcftools-1.15.1-r/bin/bcftools view --max-alleles 2 --exclude-types indels $file.vcf > $file.filt.vcf
#
#echo $file
#/programs/htslib-1.11/bin/bgzip -c $file.filt.vcf \
#    > $file.filt.vcf.gz
#/programs/htslib-1.11/bin/tabix -f -p vcf \
#    $file.filt.vcf.gz
#
