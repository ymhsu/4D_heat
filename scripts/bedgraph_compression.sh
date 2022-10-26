#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=20G
#SBATCH --mail-user=yu-ming.hsu@u-psud.fr
#SBATCH --mail-type=END
#SBATCH --job-name=bedgraph_compression
#SBATCH --array=1-19%4

source ~/.profile

feature=$(sed -n "$SLURM_ARRAY_TASK_ID"p ../data/Chip_seq/no_treatment/bw_bg/compression_list)

#compress files
gzip -c ../data/Chip_seq/no_treatment/bw_bg/sorted_dedupPCR_aboveMAPQ30_sorted_trim_M82_${feature} > ../data/Chip_seq/no_treatment/bw_bg/${feature}.gz 2> ../data/Chip_seq/no_treatment/bw_bg/${feature}.error.txt

#delete bedgraph files
rm -f ../data/Chip_seq/no_treatment/bw_bg/sorted_dedupPCR_aboveMAPQ30_sorted_trim_M82_${feature}
