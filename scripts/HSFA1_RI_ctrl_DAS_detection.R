#This work is to compare the difference between HSFA1 mutants and WT

library(tidyverse)

comp <- c("WT_HS0_vs_WT_HS1", "hsfa1a_line2_HS0_vs_hsfa1a_line2_HS1", "hsfa1a_line1_HS0_vs_hsfa1a_line1_HS1")

input_path_data_HSFA1A <-
  str_c("./data/rMATS_HSFA1A_out/", rep(comp, each = 4), "/", rep(c("RI", "SE", "A5SS", "A3SS"), length(comp)), ".MATS.JC.txt")

#create the function to import data
function_input_HSFA1A_data <- function(a){
  AS_type <- 
    str_remove(str_remove(a[1], "./data/rMATS_HSFA1A_out/.*/"), ".MATS.JC.txt")
  
  comp_HSFA1A <-
    str_remove(str_remove(a[1], "./data/rMATS_HSFA1A_out/"),
               "/.*")
  
  input_AS_data <-
    read_delim(a[1], delim = "\t", col_names = TRUE) %>%
    #https://stackoverflow.com/questions/72616662/how-to-use-stringrstr-c-with-column-indices
    mutate(ID_final = str_c(do.call(str_c, c(.[c(2, 6:11)], sep = "_")))) %>%
    mutate(AS_type = AS_type, comp = comp_HSFA1A) %>%
    mutate(AS_type = if_else(AS_type == "A5SS" | AS_type == "A3SS", str_sub(AS_type, 1, 3), AS_type)) %>%
    filter(!duplicated(ID_final))
  
  colnames(input_AS_data)[6:11] <- c("pos_1", "pos_2", "pos_3", "pos_4", "pos_5", "pos_6")
  
  input_AS_data
  
}


input_DAS_HSFA1A_combined_raw <-
  map(input_path_data_HSFA1A, function_input_HSFA1A_data) %>%
  bind_rows()

#make a table of DAS events using only padj and show their presence and absence among three lines
DAS_HSFA1A_presence_and_absence_FDR001 <-
  input_DAS_HSFA1A_combined_raw %>%
  filter(FDR < 0.01) %>%
  mutate(having_DAS = "yes") %>%
  select(ID_final, GeneID, chr, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, strand, comp, having_DAS) %>%
  pivot_wider(names_from = comp, values_from = having_DAS) %>%
  replace_na(list(WT_HS0_vs_WT_HS1 = "no", hsfa1a_line2_HS0_vs_hsfa1a_line2_HS1 = "no",
                  hsfa1a_line1_HS0_vs_hsfa1a_line1_HS1 = "no"))

write_delim(DAS_HSFA1A_presence_and_absence_FDR001, "./data/rMATS_HSFA1A_out/DAS_HSFA1A_presence_and_absence_FDR001.txt",
            delim = "\t", col_names = TRUE)

#one is available to check presence and absence of AS events among three lines
DAS_event_presence_and_absence_3_lines_statistics <-
  DAS_HSFA1A_presence_and_absence_FDR001 %>%
  group_by(WT_HS0_vs_WT_HS1, hsfa1a_line1_HS0_vs_hsfa1a_line1_HS1, hsfa1a_line2_HS0_vs_hsfa1a_line2_HS1) %>%
  summarise(count = n())

write_delim(DAS_event_presence_and_absence_3_lines_statistics, "./data/rMATS_HSFA1A_out/DAS_event_presence_and_absence_3_lines_statistics.txt",
            delim = "\t", col_names = TRUE)

### stringent filtering
min_junct=4  ## can be adjusted
padj=0.01
dPSI

input_DAS_HSFA1A_combined_stringent <-
  input_DAS_HSFA1A_combined_raw %>%
  filter(FDR < 0.01) %>%
  separate(SJC_SAMPLE_1, into=c("SJC_S1A","SJC_S1B","SJC_S1C"), sep=",")%>% ## this line may need to be adjusted depending on input data 
  separate(SJC_SAMPLE_2, into=c("SJC_S2A","SJC_S2B","SJC_S2C"), sep=",") %>%
  pivot_longer(cols = c("SJC_S1A","SJC_S1B","SJC_S1C", "SJC_S2A","SJC_S2B","SJC_S2C"),
               names_to='SJC_sample', values_to='SJC_count') %>%
  group_by(comp, AS_type, ID_final) %>%
  filter(as.double(SJC_count) > 4) %>%
  mutate(SJC_sample_count = n()) %>%
  filter(SJC_sample_count == 6) %>%
  dplyr::select(-SJC_sample, -SJC_sample_count, -SJC_count) %>%
  distinct()

input_DAS_HSFA1A_combined_stringent %>%
  group_by(comp, AS_type) %>%
  summarise(count = n())

write_delim(input_DAS_HSFA1A_combined_stringent, "./data/rMATS_HSFA1A_out/input_DAS_HSFA1A_combined_stringent_raw.txt",
            delim = "\t", col_names = TRUE)


#create TPM of all genes in different reps of HSFA1 experiments
#for grouping three lines with the similar gene expression levels

#import annotation data of all genes
M82_rMATs_anno_all_gene <- 
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name")) %>%
  mutate(gene_size = end - str) %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  filter(source == "Liftoff") %>%
  dplyr::select(seqnames = chr, start = str, end, type = feature, source, gene_name, strand) 

M82_rMATs_anno_all_gene_size <-
  M82_rMATs_anno_all_gene %>%
  mutate(size = end - start) %>%
  dplyr::select(gene_name, size) %>%
  distinct()

#import files of genes with read counts
path_read_count_gene <- 
  str_c("./data/HSFA1A_count_of_read_gene_raw/", list.files("./data/HSFA1A_count_of_read_gene_raw/", recursive = F, full.names = F))

name_HSFA1A_trt_rep <-
  str_remove(str_remove(path_read_count_gene, "./data/HSFA1A_count_of_read_gene_raw/"), "ReadsPerGene.out.tab")

read_delim(path_read_count_gene[[1]], 
           delim = "\t", col_names = c("gene_ID", "unstranded_count", "strand_count_1", "strand_count_2")) %>%
  dplyr::select(gene_ID, unstranded_count)

function_input_HSFA1A_read_count_gene <-
  function(a, b){
  input_read_count_data <-
      read_delim(a[1], 
             delim = "\t", col_names = c("gene_ID", "unstranded_count", "strand_count_1", "strand_count_2")) %>%
    dplyr::select(gene_ID, unstranded_count)
  
  colnames(input_read_count_data) <- c("gene_ID", b[1])
  
  input_read_count_data
}  

input_HSFA1A_read_count_gene_combined <-
  pmap(list(path_read_count_gene, name_HSFA1A_trt_rep), function_input_HSFA1A_read_count_gene) %>%
  reduce(., full_join)

#calculate TPM of different rep of lines among treatments
input_HSFA1A_read_count_gene_combined_TPM <-
  input_HSFA1A_read_count_gene_combined %>%
  filter(str_detect(gene_ID, "gene:")) %>%
  pivot_longer(cols = c(2:19), names_to = "trt_rep", values_to = "read_count") %>%
  replace_na(list(read_count = 0)) %>%
  #pivot_wider(names_from = "trt_rep", values_from = "read_count")
  mutate(gene_name = str_remove(gene_ID, "gene:")) %>%
  group_by(gene_name) %>%
  mutate(sum_read_count = sum(read_count)) %>%
  filter(sum_read_count != 0) %>%
  #filter(gene_name == "Solyc01g038233.1")
  mutate(label = c(1:n())) %>%
  left_join(M82_rMATs_anno_all_gene_size) %>%
  drop_na() %>%
  mutate(RPK = read_count/size) %>%
  group_by(trt_rep) %>%
  mutate(total_RPK_group = sum(RPK)/10^6) %>%
  mutate(TPM = RPK/total_RPK_group) %>%
  dplyr::select(gene_name, trt_rep, TPM) %>%
  ungroup() %>%
  mutate(trt_rep = str_c(str_replace_all(trt_rep, "-", "_"), "_TPM")) %>%
  pivot_wider(names_from = "trt_rep", values_from = "TPM")

write_delim(input_HSFA1A_read_count_gene_combined_TPM, "./data/rMATS_HSFA1A_out/input_HSFA1A_read_count_gene_combined_TPM.txt",
            delim = "\t", col_names = TRUE)

#calculate mean TPM of different rep of lines among treatments
input_HSFA1A_read_count_gene_combined_mean_TPM <-
  input_HSFA1A_read_count_gene_combined_TPM %>%
  pivot_longer(names_to = "trt_rep", values_to = "TPM", cols = c(2:19)) %>%
  mutate(trt = str_replace(trt_rep, "(.*h)_.*", "\\1")) %>%
  dplyr::select(gene_name, trt, TPM) %>%
  group_by(trt, gene_name) %>%
  summarise(mean_TPM = mean(TPM)) %>%
  ungroup() %>%
  mutate(trt = str_c(trt, "_mean_TPM")) %>%
  pivot_wider(names_from = "trt", values_from = "mean_TPM")

write_delim(input_HSFA1A_read_count_gene_combined_mean_TPM, "./data/rMATS_HSFA1A_out/input_HSFA1A_read_count_gene_combined_mean_TPM.txt",
            delim = "\t", col_names = TRUE)



#import all AS events of HSFA1A

#produce the function to import data, tidy data, and select necessary columns for AS events from 3 lines with 2 temperature conditions
produce_AS_HSFA1A <- function(a) {
  #the colnames used for later importing AS-event tables
  col_names_uni <-
    c(
      "ID_1",
      "GeneID",
      "geneSymbol",
      "chr",
      "strand",
      "pos_1",
      "pos_2",
      "pos_3",
      "pos_4",
      "pos_5",
      "pos_6",
      "ID_2",
      "IJC_SAMPLE_1",
      "SJC_SAMPLE_1",
      "IJC_SAMPLE_2",
      "SJC_SAMPLE_2",
      "IncFormLen",
      "SkipFormLen",
      "PValue",
      "FDR",
      "IncLevel1",
      "IncLevel2",
      "IncLevelDifference"
    )
  
  #set the path of AS imported from 6 directory (3 lines * 2 temperature trt)
  HSFA1A_lines <- c("hsfa1a_line1", "hsfa1a_line2", "WT")
  
  path_input_AS_raw <-
    str_c("./data/rMATS_HSFA1A_out/", rep(HSFA1A_lines, each = 2), "_HS", rep(c(0, 1), 3), "/", a[1], ".MATS.JC.txt")
  
  #import AS events
  AS_df_trt_raw <-
    map(
      path_input_AS_raw,
      read_delim,
      col_names = col_names_uni,
      delim = "\t",
      skip = 1
    )
  
  #check how many rows in each of imported AS events table
  nrow_AS_df_trt <-
    AS_df_trt_raw %>%
    map(. %>% nrow()) %>%
    as_vector()
  
  #create the combined AS tables with necessary info (chr, strand, pos, AS_type, trt, GeneID)
  AS_df_trt_raw %>%
    bind_rows() %>%
    dplyr::select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, GeneID, IJC_SAMPLE_1, SJC_SAMPLE_1, IncLevel1) %>%
    mutate(gene_name = str_remove(GeneID, "gene:")) %>%
    mutate(AS_type = a[1]) %>%
    mutate(trt = c(
      rep(path_input_AS_raw[[1]], nrow_AS_df_trt[[1]]),
      rep(path_input_AS_raw[[2]], nrow_AS_df_trt[[2]]),
      rep(path_input_AS_raw[[3]], nrow_AS_df_trt[[3]]),
      rep(path_input_AS_raw[[4]], nrow_AS_df_trt[[4]]),
      rep(path_input_AS_raw[[5]], nrow_AS_df_trt[[5]]),
      rep(path_input_AS_raw[[6]], nrow_AS_df_trt[[6]])
    )) %>%
    mutate(trt = str_remove(str_remove(trt, "./data/rMATS_HSFA1A_out/"), "/.*\\.txt")) %>%
    mutate(genotype = str_replace(trt, "(.*)_H.*", "\\1"),
           trt = str_replace(trt, ".*_(H.)*", "\\1"))
}

#import data
AS_type_import <- c("RI", "SE", "A5SS", "A3SS")


all_AS_HSFA1A_df_trt_separated_raw <-
  map(AS_type_import, produce_AS_HSFA1A) %>%
  bind_rows()

all_AS_HSFA1A_df_trt_separated_segments <- 
all_AS_HSFA1A_df_trt_separated_raw %>%
  separate(IncLevel1, into = c("PSI_s1", "PSI_s2", "PSI_s3"), sep = ",", convert = TRUE) %>%
  replace_na(list(PSI_s1 = 0, PSI_s2 = 0, PSI_s3 = 0)) %>%
  mutate(mean_PSI = c(PSI_s1 + PSI_s2 + PSI_s3)/3) %>%
  select(-IJC_SAMPLE_1, -SJC_SAMPLE_1, -PSI_s1, -PSI_s2, -PSI_s3, -GeneID) %>%
  mutate(AS_type = if_else(AS_type == "A5SS", "A5S",
                           if_else(AS_type == "A3SS", "A3S", AS_type))) %>%
  mutate(str = if_else(AS_type == "A5S" & strand == "+", pos_4, 
                       if_else(AS_type == "A5S" & strand == "-", pos_1,
                               if_else(AS_type == "A3S" & strand == "+", pos_1,
                                       if_else(AS_type == "A3S" & strand == "-", pos_4, 
                                               if_else(AS_type == "SE", pos_1, pos_4)))))) %>%
  mutate(end = if_else(AS_type == "A5S" & strand == "+", pos_2, 
                       if_else(AS_type == "A5S" & strand == "-", pos_3,
                               if_else(AS_type == "A3S" & strand == "+", pos_3,
                                       if_else(AS_type == "A3S" & strand == "-", pos_2, 
                                               if_else(AS_type == "SE", pos_2, pos_5)))))) %>%
  dplyr::select(seqnames = chr, start = str, end, AS_type, gene_name, strand, genotype, trt, mean_PSI) %>%
  arrange(AS_type, trt, seqnames, start) 


write_delim(all_AS_HSFA1A_df_trt_separated_segments, "./data/rMATS_HSFA1A_out/all_AS_HSFA1A_df_trt_separated_segments.bed", delim = "\t", col_names = TRUE)

input_DAS_HSFA1A_combined_stringent <-
  read_delim("./data/rMATS_HSFA1A_out/input_DAS_HSFA1A_combined_stringent_raw.txt",
            delim = "\t", col_names = TRUE)

input_DAS_HSFA1A_combined_stringent_bed6 <- 
  input_DAS_HSFA1A_combined_stringent %>%
  mutate(AS_type = if_else(AS_type == "A5SS", "A5S",
                           if_else(AS_type == "A3SS", "A3S", AS_type))) %>%
  mutate(str = if_else(AS_type == "A5S" & strand == "+", pos_4, 
                       if_else(AS_type == "A5S" & strand == "-", pos_1,
                               if_else(AS_type == "A3S" & strand == "+", pos_1,
                                       if_else(AS_type == "A3S" & strand == "-", pos_4, 
                                               if_else(AS_type == "SE", pos_1, pos_4)))))) %>%
  mutate(end = if_else(AS_type == "A5S" & strand == "+", pos_2, 
                       if_else(AS_type == "A5S" & strand == "-", pos_3,
                               if_else(AS_type == "A3S" & strand == "+", pos_3,
                                       if_else(AS_type == "A3S" & strand == "-", pos_2, 
                                               if_else(AS_type == "SE", pos_2, pos_5)))))) %>%
  mutate(gene_name = str_remove(GeneID, "gene:")) %>%
  dplyr::select(seqnames = chr, start = str, end, AS_type, gene_name, strand, comp) %>%
  mutate(AS_group_type = "DAS", genotype = str_replace(comp, "(.*)_HS0_vs.*", "\\1")) %>%
  arrange(AS_type, seqnames, start)

input_DAS_HSFA1A_combined_stringent_bed6 %>%
  View()
all_HSFA1A_AS_DAS_combined_table_raw <-
all_AS_HSFA1A_df_trt_separated_segments %>%
  left_join(input_DAS_HSFA1A_combined_stringent_bed6) %>%
  replace_na(list(comp = "no_comp", AS_group_type = "non_DAS_AS"))


all_HSFA1A_AS_DAS_combined_table_raw %>%
  group_by(genotype, trt, AS_type) %>%
  summarise(count = n()) %>%
  View()


#produce genes containing AS (stable PSI or DAS)
HSFA1A_all_genes_overlapped_AS <- 
  all_HSFA1A_AS_DAS_combined_table_raw %>%
  mutate(overlap_AS = "yes") %>%
  dplyr::select(gene_name, genotype, overlap_AS) %>%
  distinct()

HSFA1A_all_genes_overlapped_AS %>%
  group_by(genotype) %>%
  summarise()

#import all introns from annotation data to produce the set of ctrl introns 
#import all introns
M82_rMATs_anno_all_intron_exon <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron_exon_ordered_corrected.bed", delim = "\t") %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  dplyr::select(-genetic_feature_label) %>%
  filter(source != "REGARN")

M82_rMATs_anno_all_intron <-
  M82_rMATs_anno_all_intron_exon %>%
  filter(feature == "intron")

#https://stackoverflow.com/questions/19590541/r-duplicate-a-matrix-several-times-and-then-bind-by-rows-together

#produce all introns of 3 AS groups from genes larger than 3.5 kb (including all types of AS)
all_AS_DAS_ctrl_introns_HSFA1A_raw <-
  do.call(bind_rows, replicate(3, M82_rMATs_anno_all_intron, simplify=FALSE)) %>%
  mutate(genotype = rep(c("WT", "hsfa1a_line1", "hsfa1a_line2"), each = nrow(M82_rMATs_anno_all_intron))) %>%
  left_join(HSFA1A_all_genes_overlapped_AS) %>%
  replace_na(list(overlap_AS = "no")) %>%
  filter(overlap_AS != "yes") %>%
  mutate(AS_type = "noAS", trt = "all", comp = "no_comp", AS_group_type = "ctrl", mean_PSI = 0) %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, genotype, trt, mean_PSI, comp, AS_group_type) %>%
  bind_rows(., all_HSFA1A_AS_DAS_combined_table_raw) %>%
  mutate(chr = str_remove(seqnames, "chr")) %>%
  filter(str_detect(chr, "[a-z]", negate = TRUE)) %>%
  mutate(chr = as.double(chr)) %>%
  drop_na() %>%
  arrange(AS_group_type, AS_type, chr, start) %>%
  left_join(M82_rMATs_anno_all_gene_size) %>%
  filter(size > 3500) %>%
  dplyr::select(-size) 

write_delim(all_AS_DAS_ctrl_introns_HSFA1A_raw, "./data/rMATS_HSFA1A_out/all_AS_DAS_ctrl_introns_HSFA1A_raw.bed", delim = "\t", col_names = TRUE)

#produce all introns of 3 AS groups from genes larger than 3.5 kb (including only RI, the comp between HS1 and HS6 is removed from DAS)
all_RI_stable_PSI_DAS_ctrl_introns_HSFA1A_raw <-
  all_AS_DAS_ctrl_introns_HSFA1A_raw %>%
  filter(AS_type == "RI" | AS_type == "noAS") %>%
  #filter(comp != "HS1_HS6") %>%
  dplyr::select(-trt, -mean_PSI, -comp) %>%
  distinct() 

write_delim(all_RI_stable_PSI_DAS_ctrl_introns_HSFA1A_raw, "./data/rMATS_HSFA1A_out/all_RI_stable_PSI_DAS_ctrl_introns_splicing_sites_HSFA1A_raw.bed", delim = "\t", col_names = TRUE)


all_RI_stable_PSI_DAS_ctrl_introns_HSFA1A_raw %>%
  group_by(genotype, AS_group_type) %>%
  summarise(count = n())

#import TPM of three lines and calculaute quantiles based on their HS0 TPM
input_HSFA1A_read_count_gene_combined_mean_TPM <-
  read_delim("./data/rMATS_HSFA1A_out/input_HSFA1A_read_count_gene_combined_mean_TPM.txt",
            delim = "\t", col_names = TRUE)

all_HSFA1A_gene_expression_level_0_removed <-
  input_HSFA1A_read_count_gene_combined_mean_TPM %>%
  pivot_longer(cols = c(2:7), names_to = "name_TPM_trt", values_to = "mean_TPM") %>%
  filter(mean_TPM != 0) %>%
  group_by(name_TPM_trt) %>%
  mutate(quantile_group = if_else(mean_TPM < quantile(mean_TPM, 0.1), "q1",
                                  if_else(quantile(mean_TPM, 0.1) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.2), "q2",
                                          if_else(quantile(mean_TPM, 0.2) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.3), "q3",
                                                  if_else(quantile(mean_TPM, 0.3) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.4), "q4",
                                                          if_else(quantile(mean_TPM, 0.4) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.5), "q5",
                                                                  if_else(quantile(mean_TPM, 0.5) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.6), "q6",
                                                                          if_else(quantile(mean_TPM, 0.6) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.7), "q7",
                                                                                  if_else(quantile(mean_TPM, 0.7) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.8), "q8",
                                                                                          if_else(quantile(mean_TPM, 0.8) <= mean_TPM & mean_TPM < quantile(mean_TPM, 0.9), "q9", "q10")))))))))) 
#extract quantile group of 0H of 3 lines  
all_HSFA1A_gene_expression_level_0_removed_0H <-
  all_HSFA1A_gene_expression_level_0_removed %>%
  mutate(genotype = str_replace(name_TPM_trt, "MM_(.*)_HS.*", "\\1")) %>%
  filter(str_detect(name_TPM_trt, "0h")) %>%
  ungroup() %>%
  dplyr::select(gene_name, genotype, quantile_group) %>%
  distinct()

RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw <-
  all_RI_stable_PSI_DAS_ctrl_introns_HSFA1A_raw %>%
  left_join(all_HSFA1A_gene_expression_level_0_removed_0H) %>%
  drop_na() %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, genotype, AS_group_type, chr, quantile_group) %>%
  arrange(genotype, chr, start)

write_delim(RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw, "./data/rMATS_HSFA1A_out/RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw.bed", delim = "\t", col_names = TRUE)


RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw %>%
  group_by(genotype, AS_group_type, quantile_group) %>%
  summarise(count = n()) %>%
  View()

genotype_HSFA1A <-
  unique(RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw$genotype)

#produce files for quantiles by event numbers (splicing sites)
produce_norm_expression_level_RI_AS_DAS_HSFA1 <-
  function(a, b){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw %>%
      filter(quantile_group == a[1] & genotype == b[1]) %>%
      split(.$AS_group_type)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level_HSFA1/", b[1], "_", names(output_l), "_", output_l[[1]]$quantile_group[[1]], "_ss_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

pwalk(list(rep(str_c("q", c(1:10)), 3), rep(genotype_HSFA1A, each = 10)), produce_norm_expression_level_RI_AS_DAS_HSFA1)

#produce files for quantiles by event numbers (gene bodies)
RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_gene_bodies_raw <-
  RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_raw %>%
  dplyr::select(AS_type, gene_name, strand, genotype, AS_group_type, chr, quantile_group) %>%
  left_join(M82_rMATs_anno_all_gene) %>%
  distinct() %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, genotype, AS_group_type, chr, quantile_group) %>%
  arrange(genotype, chr, start) %>%
  group_by(genotype, AS_group_type, quantile_group) %>%
  arrange(genotype, chr, start)

produce_norm_expression_level_RI_AS_DAS_HSFA1_gene_bodies <-
  function(a, b){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_HSFA1A_gene_bodies_raw %>%
      filter(quantile_group == a[1] & genotype == b[1]) %>%
      split(.$AS_group_type)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level_HSFA1/", b[1], "_", names(output_l), "_", output_l[[1]]$quantile_group[[1]], "_gene_bodies_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

pwalk(list(rep(str_c("q", c(1:10)), 3), rep(genotype_HSFA1A, each = 10)), produce_norm_expression_level_RI_AS_DAS_HSFA1_gene_bodies)


