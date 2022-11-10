#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --mail-user=yu-ming.hsu@u-psud.fr
#SBATCH --mail-type=END
#SBATCH --job-name=HS0_bamcompare_AS_blacklist
#SBATCH --array=1-19

source ~/.profile
module load deeptools

marks_name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../data/Chip_seq/no_treatment/bam/mark_list)
marks_f=/NetScratch/REGARN/yhsu/4D_heat/data/Chip_seq/no_treatment/bam/sorted_dedupPCR_aboveMAPQ30_sorted_trim_M82_${marks_name}.bam
ctrl=/NetScratch/REGARN/yhsu/4D_heat/data/Chip_seq/no_treatment/bam/sorted_dedupPCR_aboveMAPQ30_sorted_trim_M82_Input_S12_R1_001_bowtie_res.bam


for i in $(cat ../data/AS_bed_for_ML/AS_blacklist_list)
do

AS_label=$(echo ${i} | sed 's/_for_ML//g' | sed 's/_blacklist.bed//g')

#echo ${test} >> output_test
bamCompare -b1 ${marks_f} \
-b2 ${ctrl} \
-o /NetScratch/REGARN/yhsu/4D_heat/data/AS_bed_for_ML/${marks_name}_${AS_label}.bedgraph \
-of bedgraph \
--blackListFileName ../data/AS_bed_for_ML/${i} \
--operation ratio \
--binSize 1 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 5 2> ../data/AS_bed_for_ML/${marks_name}_${AS_label}_bamCompare.log

done

module unload deeptools 

