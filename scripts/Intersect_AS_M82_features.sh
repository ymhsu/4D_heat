#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=10G
#SBATCH --mail-user=yu-ming.hsu@u-psud.fr
#SBATCH --mail-type=END
#SBATCH --job-name=Intersect_AS_M82_features


source ~/.profile
module load bedtools2

ctrl="AS_control AS_sig_FDR05"
comp="HS6_HS0 HS1_HS0 HS1_HS6"
AS="A5S A3S RI SE"
feature="gene exon_order intron_order"


for s in ${ctrl}
do
for i in ${comp}
do
for j in ${AS}
do
for k in ${feature}
do

bedtools intersect -a ../data/AS_control_set/${s}_TPM_q05_${i}_${j}.bed \
-b ../data/M82_annotation_data/M82_rMATs_anno_all_${k}.bed -wb \
| bedtools sort > ../data/Intersected_AS_TPM_q05_M82_anno_rMATs/${s}_TPM_q05_${i}_${j}_intersected_M82_rMATs_${k}.bed \
2> ../data/Intersected_AS_TPM_q05_M82_anno_rMATs/${s}_TPM_q05_${i}_${j}_intersected_M82_rMATs_${k}.error

done
done
done
done

module unload bedtools2
