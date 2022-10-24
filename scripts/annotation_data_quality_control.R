library(tidyverse)

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
