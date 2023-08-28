library(tidyverse)

#extract DEG (non DAS) vs DAS (non-DE) to compare their peak profiles

read_DEseq_data <- function(a){
  read_delim(a[1], delim = "\t", skip = 1,
             col_names = c("gene_name", "baseMean", "log2FoldChange", "lfcSE", "stat",
                           "pvalue", "padj"))
}

path_to_DEseq_output <- str_c("./data/DEseq_result/0h", c(1, 6), "h_all_data.tsv")

#extract DEGs are defined as genes with a change in expression of log2 Fold Change |1| and adjusted p value below 0.05
#https://www-ncbi-nlm-nih-gov.insb.bib.cnrs.fr/pmc/articles/PMC5452352/
DEseq_output_raw <-
  map(path_to_DEseq_output, read_DEseq_data) %>%
  map(. %>% filter(pvalue < 0.05 & abs(log2FoldChange) >= 1)) %>%
  map(. %>% mutate(gene_name = str_remove(gene_name, "gene:")))

DEG_gene_list_0_1 <-
  DEseq_output_raw[[1]] %>%
  dplyr::select(gene_name) %>%
  mutate(DEG_type = "DEG_0_1")

DEG_gene_list_0_6 <-
  DEseq_output_raw[[2]] %>%
  dplyr::select(gene_name) %>%
  mutate(DEG_type = "DEG_0_6")

#import the file of refined DAS genes and remove DAS genes from DEG
all_AS_DAS_ctrl_introns_raw <-
  read_delim("./data/AS_df_trt/all_AS_DAS_ctrl_introns_raw.bed", delim = "\t", col_names = TRUE)

DAS_gene_list_all_heat_comp <-
  all_AS_DAS_ctrl_introns_raw %>%
  filter(group_type == "DAS") %>%
  dplyr::select(gene_name, AS_group_type = group_type) %>%
  distinct()

DEG_gene_non_DAS_01_06 <-
  bind_rows(DEG_gene_list_0_1, DEG_gene_list_0_6) %>%
  left_join(DAS_gene_list_all_heat_comp) %>%
  replace_na(list(AS_group_type = "non_DAS")) %>%
  filter(AS_group_type == "non_DAS") %>%
  group_by(gene_name) %>%
  mutate(comp = "no_heat_comp")

#extract DEG from DAS gene
DEG_gene_01_06 <-
  bind_rows(DEG_gene_list_0_1, DEG_gene_list_0_6)

all_RI_DAS_non_DEG_gene_list <-
  all_AS_DAS_ctrl_introns_raw %>%
  filter(group_type == "DAS") %>%
  left_join(DEG_gene_01_06) %>%
  replace_na(list(DEG_type = "non_DEG")) %>%
  filter(str_detect(DEG_type, "non_DEG")) %>%
  distinct() %>%
  filter(AS_type == "RI") %>%
  dplyr::select(gene_name, DEG_type, AS_group_type = group_type, comp) %>%
  distinct() %>%
  group_by(comp) 

#combined the information of DEG (non-DAS) and DAS (non-DEG) with the coordinate from annotation data
#and also the expression level at HS0 
#import the information of gene size for normalizing the peak signal (similar to RPKM)
M82_rMATs_anno_all_gene <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name")) %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  filter(source == "Liftoff")

#import RNA expression level files to compare expression levels of no-AS, AS (stable PSI) and DAS groups
expression_level_raw <- read_csv("./data/expression_level.csv", col_names = TRUE) %>%
  mutate(ID = str_remove(ID, "gene:"))

expression_level <- expression_level_raw %>%
  mutate(mean_0H = (RNA_M82_HS_0h_rep1_S4 + RNA_M82_HS_0h_rep2_S5)/2,
         mean_1H = (RNA_M82_HS_1h_rep1_S6 + RNA_M82_HS_1h_rep2_S7)/2,
         mean_6H = (RNA_M82_HS_6h_rep1_S8 + RNA_M82_HS_6h_rep2_S9)/2) %>%
  dplyr::select(gene_name = ID, mean_0H, mean_1H, mean_6H)

expression_level_with_quantile <-
M82_rMATs_anno_all_gene %>%
  left_join(expression_level) %>%
  drop_na() %>%
  filter(mean_0H != 0) %>%
  mutate(quantile_group_0H = if_else(mean_0H < quantile(mean_0H, 0.1), "q1",
                                     if_else(quantile(mean_0H, 0.1) <= mean_0H & mean_0H < quantile(mean_0H, 0.2), "q2",
                                             if_else(quantile(mean_0H, 0.2) <= mean_0H & mean_0H < quantile(mean_0H, 0.3), "q3",
                                                     if_else(quantile(mean_0H, 0.3) <= mean_0H & mean_0H < quantile(mean_0H, 0.4), "q4",
                                                             if_else(quantile(mean_0H, 0.4) <= mean_0H & mean_0H < quantile(mean_0H, 0.5), "q5",
                                                                     if_else(quantile(mean_0H, 0.5) <= mean_0H & mean_0H < quantile(mean_0H, 0.6), "q6",
                                                                             if_else(quantile(mean_0H, 0.6) <= mean_0H & mean_0H < quantile(mean_0H, 0.7), "q7",
                                                                                     if_else(quantile(mean_0H, 0.7) <= mean_0H & mean_0H < quantile(mean_0H, 0.8), "q8",
                                                                                             if_else(quantile(mean_0H, 0.8) <= mean_0H & mean_0H < quantile(mean_0H, 0.9), "q9", "q10")))))))))) %>%
  dplyr::select(gene_name, mean_0H, mean_1H, mean_6H, quantile_group_0H)

#check the gene count of non-DAS DEG and non-DEG DAS according to different quantiles defined by expression levels of HS0
combined_non_DEG_DAS_and_non_DAS_DEG_list <-
  bind_rows(DEG_gene_non_DAS_01_06, all_RI_DAS_non_DEG_gene_list) %>%
  left_join(expression_level_with_quantile) %>%
  filter(quantile_group_0H != "q1" & quantile_group_0H != "q2" & quantile_group_0H != "q3")

combined_non_DEG_DAS_and_non_DAS_DEG_list %>%
  group_by(DEG_type, AS_group_type, comp, quantile_group_0H) %>%
  summarise(mean_ex_l_HS0 = median(mean_0H), mean_ex_l_HS1 = median(mean_1H), count = n()) %>%
  View()

raw_output_non_DEG_DAS_and_non_DAS_DEG_list <-
  M82_rMATs_anno_all_gene %>%
  left_join(combined_non_DEG_DAS_and_non_DAS_DEG_list) %>%
  drop_na() %>%
  filter(end - str > 3500) %>%
  filter(comp != "HS1_HS6") %>%
  dplyr::select(chr, str, end, gene_name, feature, strand, DEG_type, AS_group_type, comp, quantile_group_0H)

raw_output_non_DEG_DAS_and_non_DAS_DEG_list %>%
  group_by(DEG_type, AS_group_type, comp, quantile_group_0H) %>%
  #filter(end - str > 3500) %>%
  summarise(median_gene_size = median(end - str), count = n()) %>%
  View()
  
  
"DEG_0_1" %in% raw_output_non_DEG_DAS_and_non_DAS_DEG_list$DEG_type

create_non_DEG_DAS_and_non_DAS_DEG_list_with_q_ex_l <- function(a){
if (a[1] %in% raw_output_non_DEG_DAS_and_non_DAS_DEG_list$DEG_type) {
  non_DAS_DEG_list <-
  raw_output_non_DEG_DAS_and_non_DAS_DEG_list %>%
    filter(DEG_type == a[1]) %>%
    split(.$quantile_group_0H)
  
  output_quantile_names <- str_c("./data/gene_bodies_non_DEG_DAS_non_DAS_DEG/gene_bodies_non_DAS_", a[1], "_ex_l_", names(non_DAS_DEG_list), "_bed6.bed")
  
  pwalk(list(non_DAS_DEG_list, output_quantile_names), write_delim, delim = "\t", col_names = FALSE)
  
} else {
  non_DEG_DAS_list <-
  raw_output_non_DEG_DAS_and_non_DAS_DEG_list %>%
    filter(comp == a[1])%>%
    split(.$quantile_group_0H)
  
  output_quantile_names <- str_c("./data/gene_bodies_non_DEG_DAS_non_DAS_DEG/gene_bodies_non_DEG_DAS_", a[1], "_ex_l_", names(non_DEG_DAS_list), "_bed6.bed")
  
  pwalk(list(non_DEG_DAS_list, output_quantile_names), write_delim, delim = "\t", col_names = FALSE)
}
}

walk(c("DEG_0_1", "DEG_0_6", "HS0_HS1", "HS0_HS6"), create_non_DEG_DAS_and_non_DAS_DEG_list_with_q_ex_l)  

