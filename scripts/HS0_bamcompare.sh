#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --mail-user=yu-ming.hsu@u-psud.fr
#SBATCH --mail-type=END
#SBATCH --job-name=test_bamcompare
#SBATCH --array=1-19

source ~/.profile
module load deeptools samtools

marks_name=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../data/Chip_seq/no_treatment/bam/mark_list)
marks_f=/NetScratch/REGARN/yhsu/4D_heat/data/Chip_seq/no_treatment/bam/${marks_name}.bam
ctrl=/NetScratch/REGARN/yhsu/4D_heat/data/Chip_seq/no_treatment/bam/sorted_dedupPCR_aboveMAPQ30_sorted_trim_M82_Input_S12_R1_001_bowtie_res.bam


samtools index ${marks_f}

bamCompare -b1 ${marks_f} \
-b2 ${ctrl} \
-o /NetScratch/REGARN/yhsu/4D_heat/data/Chip_seq/no_treatment/bw_bg/${marks_name}.bedgraph \
-of bedgraph \
--binSize 1 \
--smoothLength 60 \
--extendReads 150 \
--centerReads \
-p 5 2> ../data/Chip_seq/no_treatment/bw_bg/${marks_name}_bamCompare.log


module unload deeptools samtools
