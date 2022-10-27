library(tidyverse)

#this is to check some info related to gene, exons and mRNAs
M82_gene <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12Chr_gene.bed", delim = "\t", col_names = c("Chr", "str", "end", "info"))
M82_mRNA <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12Chr_mRNA.bed", delim = "\t", col_names = c("Chr", "str", "end", "info"))
M82_exons <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12Chr_exon.bed", delim = "\t", col_names = c("Chr", "str", "end", "info"))

?str_detect
location_mRNA_exon <- str_locate(M82_exons$info, "Parent=mRNA:")
location_copy_exon <- str_locate(M82_exons$info, ";extra_copy")
str(location_mRNA_exon)

location_mRNA_exon[,2]
M82_exons$info

M82_exons_test <- M82_exons %>%
  mutate(ext_str = location_mRNA_exon[,2] + 1, ext_end = location_copy_exon[,1] - 1) %>%
  mutate(mRNA = str_sub(info, ext_str, ext_end)) %>%
  select(Chr, str, end, mRNA) 

M82_exons_test_light <- M82_exons_test %>%
  group_by(Chr, mRNA) %>%
  summarise(str = min(str), end = max(end)) %>%
  arrange(Chr, str, end, mRNA) %>%
  ungroup() %>%
  select(Chr, str, end)

M82_mRNA_light <- M82_mRNA %>%
  select(Chr, str, end) %>%
  arrange(Chr, str, end)

M82_gene_light <- M82_gene %>%
  select(Chr, str, end) %>%
  arrange(Chr, str, end)

setdiff(M82_gene_light, M82_mRNA_light)
intersect(M82_gene_light, M82_mRNA_light)

#define introns using the annotation data defined by Jeremie in rMATs

M82_rMATs_anno_all <- read_delim("./data/M82_annotation_data/Sl_M82_20220504_all.gff3", skip = 3, col_names = c("chr", "source", "feature", "str", "end", "no_1", "strand", "cod", "anno")) %>%
  select(chr, str, end, strand, feature, source, anno)

loc_M82_rMATs_anno_end <- str_locate(M82_rMATs_anno_all$anno, ";gene_id") 
loc_M82_rMATs_anno_end[,1]

?str_remove

M82_rMATs_anno_all %>%
  filter(source == "REGARN") %>%
  select(chr) %>%
  distinct()

#produce the bed file of gene using the annotation data of Jeremie
M82_rMATs_anno_all_exon_bed_l <- M82_rMATs_anno_all %>%
  mutate(loc_M82_rMATs_anno_end = loc_M82_rMATs_anno_end[, 1]) %>%
  mutate(anno = str_sub(anno, 1, loc_M82_rMATs_anno_end-1)) %>%
  mutate(anno = str_remove(anno, "transcript_id.*:")) %>%
  mutate(anno = str_remove(anno, "transcript_id=")) %>%
  select(-loc_M82_rMATs_anno_end) %>%
  filter(feature == "exon") %>%
  filter(str_detect(chr, "chr")==TRUE) %>%
  mutate(chr_number = as.double(str_sub(chr, 4, 5)), str = str - 1) %>%
  filter(chr_number>=1 & chr_number <= 12) %>%
  arrange(chr_number, str, end) %>%
  select(-chr_number) %>%
  split(.$source)

M82_rMATs_anno_all_gene_bed_l <- M82_rMATs_anno_all_exon_bed_l %>%
  bind_rows() %>%
  mutate(chr_number = as.double(str_sub(chr, 4, 5))) %>%
  group_by(source, anno, chr, chr_number, strand) %>%
  summarise(str = min(str), end = max(end)) %>%
  mutate(feature = "gene") %>%
  select(chr, str, end, strand, feature, source, anno, chr_number) %>%
  ungroup() %>%
  arrange(chr_number, str, end) %>%
  select(-chr_number) %>%
  split(.$source)

M82_rMATs_anno_all_gene_bed_l[[1]]

write_delim(M82_rMATs_anno_all_gene_bed_l[[1]], "./data/M82_annotation_data/M82_rMATs_anno_all_gene_liftoff.bed", col_names = FALSE, delim = "\t")
write_delim(M82_rMATs_anno_all_gene_bed_l[[2]], "./data/M82_annotation_data/M82_rMATs_anno_all_gene_regarn.bed", col_names = FALSE, delim = "\t")

write_delim(M82_rMATs_anno_all_exon_bed_l[[1]], "./data/M82_annotation_data/M82_rMATs_anno_all_exon_liftoff.bed", col_names = FALSE, delim = "\t")
write_delim(M82_rMATs_anno_all_exon_bed_l[[2]], "./data/M82_annotation_data/M82_rMATs_anno_all_exon_regarn.bed", col_names = FALSE, delim = "\t")
