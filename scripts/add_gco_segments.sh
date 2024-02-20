#!/bin/bash
#SBATCH --partition long7
#SBATCH --nodes 1
#SBATCH --ntasks-per-node 1
#SBATCH --mem-per-cpu 8G
#SBATCH --job-name add_gco_segments
#SBATCH --output /home/sna53/gco_logs/add_gco_segments.%A_%a.log
#SBATCH --array=1-4400

/programs/bin/labutils/mount_server cbsubscb18 /storage
cd /fs/cbsubscb18/storage/siddharth/gene-conversions

counter=1
for chr in {1}
do
	for hap in {1..200}
	do
		if [ $counter -eq $SLURM_ARRAY_TASK_ID ]
		then
			break 2
		fi
		counter=$(($counter+1))
	done
done

GCO_DIR=v2_data/d00_inputs/gcos/${chr}

SOURCE_DIR=v2_data/d00_inputs/source_data

function roll_dice () {
    min=1
    max=$1
    number=$(expr $min + $RANDOM % $max)
    echo $number
}

#Get gco list for current individual
sed -n "${hap}p" ${GCO_DIR}/gcos_all.gco_chr${chr}.asu.mrk.bp | tr " " "\n" | tail -n +2 \
	> ${GCO_DIR}/gcos.chr${chr}.hap${hap}.list.txt
gcolist=$(cat ${GCO_DIR}/gcos.chr${chr}.hap${hap}.list.txt)

#Split admixed files into individual haplotype phgenos
awk -v i=$hap -F '' '{print $i}' ${GCO_DIR}/gcos_all.recomb_chr${chr}.phgeno \
	> ${GCO_DIR}/gcos_all.recomb_chr${chr}.hap${hap}.phgeno

#Get gco segments and splice it into admixed individual phgeno
main=${GCO_DIR}/gcos_all.recomb_chr${chr}.hap${hap}.phgeno
gco_num=1
echo "size of original phgeno = " && wc -l $main

touch ${GCO_DIR}/gcos_table_hap${hap}.txt
for line in $gcolist
do
    pop=$(awk -F"/|:" '{print $1}' <(echo $line))
	start=$(awk -F"/|:" '{print $2}' <(echo $line))
	stop=$(($(awk -F"/|:" '{print $3}' <(echo $line))-1))
	quit=$(($stop+1))
	echo "start , stop, quit = $start , $stop, $quit"
    

	if [ $pop -eq 0 ]
	then
		echo "pop is YRI"
		panelfile=${SOURCE_DIR}/OutOfAfrica_4J17_YRI_chr${chr}.sim.phgeno
	elif [ $pop -eq 1 ]
	then
		echo "pop is CEU"
		panelfile=${SOURCE_DIR}/OutOfAfrica_4J17_CEU_chr${chr}.sim.phgeno
	else
		echo "error invalid population"
	fi

	#Get gco segment
	rand_hap=$(roll_dice 240)
	echo "sampling from hap : $rand_hap"
	touch ${GCO_DIR}/gcos.chr${chr}.hap${hap}.gco${gco_num}.phgeno
	sed -n "${start},${stop}p;${quit}q" <(awk -v i=$rand_hap -F '' '{print $i}' $panelfile)\
		> ${GCO_DIR}/gcos.chr${chr}.hap${hap}.gco${gco_num}.phgeno


	#splice gco segment into the recomb phgeno file
	gco=${GCO_DIR}/gcos.chr${chr}.hap${hap}.gco${gco_num}.phgeno
	
	gco_num=$(($gco_num+1)) #increment right before writing the spliced file
	
    if [ $start -lt $stop ] #Splice file only if gco has at least one marker or else it will leave an empty line
    then
        {
	    head -n $((start-1)) $main
	    cat $gco
	    tail -n $(($(wc -l < $main)-$stop)) $main
	    } > ${GCO_DIR}/gcos_all.recomb.gco_chr${chr}.hap${hap}.gco${gco_num}.phgeno
    
    else
        echo "start $start is not strictly less than stop $stop - skipping splicing operation"
        cat $main > ${GCO_DIR}/gcos_all.recomb.gco_chr${chr}.hap${hap}.gco${gco_num}.phgeno
    fi

    #print status after splicing	
    echo "size of edited phgeno = $(wc -l \
		< ${GCO_DIR}/gcos_all.recomb.gco_chr${chr}.hap${hap}.gco${gco_num}.phgeno)"
    num_diffs=$(diff -U 0 \
		${GCO_DIR}/gcos_all.recomb.gco_chr${chr}.hap${hap}.gco${gco_num}.phgeno $main \
		| grep ^@ | wc -l )
	echo "number of differences = $num_diffs"	
	main=${GCO_DIR}/gcos_all.recomb.gco_chr${chr}.hap${hap}.gco${gco_num}.phgeno
    
    #Write to table of ground truths
    echo -e "$hap\t$chr\t$pop\t$start\t$stop\t$num_diffs" \
        >> ${GCO_DIR}/gcos_table_hap${hap}.txt
done
