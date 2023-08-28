#This script is to use refined peak of epigenetic features and RI events to run boruta
#We only perform pairwise comparisons between different PSI groups (stable-PSI and DAS are combined together)
#The refined peak will be based on q-value below 0.001
#And the events will be based on normalized by their expression levels at HS0

#import packages
Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges", "Boruta")
lapply(Packages, library, character.only = TRUE)


#import all splicing sites of 3 AS groups (after refinement of their gene sizes, larger than 3.5 kb)
#and make the list of RI events
all_RI_stable_PSI_DAS_ctrl_introns_raw <-
  read_delim("./data/AS_df_trt/all_RI_stable_PSI_DAS_ctrl_introns_splicing_sites_raw.bed", delim = "\t", col_names = TRUE) %>%
  filter(group_type != "ctrl") %>%
  dplyr::select(seqnames, start, end) %>%
  mutate(having_RI = "yes")

#import the information of gene size for normalizing the peak signal (similar to RPKM)
M82_rMATs_anno_all_gene <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name")) %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  filter(source == "Liftoff")

M82_all_gene_size <-
  M82_rMATs_anno_all_gene %>%
  dplyr::select(seqnames = chr, start = str, end, type = feature, source, gene_name, strand) %>%
  mutate(size = end - start) %>%
  dplyr::select(gene_name, size)

#import all RI events from HS0 stage
#define the name of the table imported later
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

RI_HS0_raw <- read_delim(str_c("./data/AS_df_trt/HS0/RI.MATS.JC.txt"), col_names = col_names_uni, skip = 1) %>%
  dplyr::select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, IncLevel1, GeneID) %>%
  filter(str_detect(IncLevel1, "NA")==FALSE) %>%
  mutate(AS_type = "RI") %>%
  mutate(GeneID = str_remove(GeneID, "gene:"))



#split IncLevel1 (there are two replicates which will be used for calculating average PSI)
RI_HS0_raw_IncLevel1_m <- str_split(RI_HS0_raw$IncLevel1, ",", simplify = TRUE)

#add PSI data from two replicates to compute "PSI" which is the average PSI of r1 and r2 
RI_HS0_PSI <-
  RI_HS0_raw %>%
  mutate(
    IncLevel_r1 = as.double(RI_HS0_raw_IncLevel1_m[,1]),
    IncLevel_r2 = as.double(RI_HS0_raw_IncLevel1_m[,2])
  ) %>%
  mutate(PSI = (IncLevel_r1 + IncLevel_r2)/2) %>%
  dplyr::select(-IncLevel_r1, -IncLevel_r2) %>%
  mutate(AS_type = if_else(AS_type == "A5SS", "A5S",
                           if_else(AS_type == "A3SS", "A3S", AS_type))) %>%
  mutate(str = pos_4, end = pos_5) %>%
  dplyr::select(seqnames = chr, start = str, end, AS_type, gene_name = GeneID, strand, PSI) %>%
  arrange(seqnames, start) %>%
  left_join(all_RI_stable_PSI_DAS_ctrl_introns_raw) %>%
  drop_na() %>%
  mutate(PSI_group = if_else(PSI <= 0.05, "PSI_5",
                             if_else(PSI > 0.05 & PSI <= 0.2, "PSI_5_20",
                                     if_else(PSI > 0.2 & PSI <= 0.4, "PSI_20_40", 
                                             if_else(PSI > 0.4 & PSI <= 0.6, "PSI_40_60",
                                                     if_else(PSI > 0.6 & PSI <= 0.8, "PSI_60_80", 
                                                             if_else(PSI > 0.8 & PSI <= 0.95, "PSI_80_95", "PSI_95"))))))) %>%
  left_join(M82_all_gene_size) %>%
  distinct()

write_delim(RI_HS0_PSI, "./data/refined_AS_bed_for_ML_HS0/refined_all_RI_HS0_PSI_events.txt", delim = "\t", col_names = TRUE)

#import the signal of 20 epigenetic marks for joining their information of peak signals to RI events
epimark_HS0 <- c("Pol2", "H3K14ac", "H3K4ac", "H3K18ac", "H3K27ac",
                 "H3K27me3", "H3K4me1", "H3K4me3", "H3K9ac", "ATAC",
                 "H3K36me3", "H3K79ac", "H4K12ac", "H4K16ac", "H4K20ac",
                 "H4K5ac", "H4K8ac", "H3K36ac", "H3K4me2", "H3K56ac")

path_epimark_q001 <- str_c("./data/peak_epimarks/combined_epimark_fold_enrichment_q001/M82_0h_", epimark_HS0,
                           "_fold_enrichment_q001.bedgraph")
import_epimark <-
  function(a){
  read_delim(str_c("./data/peak_epimarks/combined_epimark_fold_enrichment_q001/M82_0h_", a[1],
                   "_fold_enrichment_q001.bedgraph"), delim = "\t", col_names = c("seqnames", "start", "end", "fold_enrichment", "name")) %>%
    mutate(epimark = a[1])
}

refined_epimark <-
  map(epimark_HS0, import_epimark)

#since we use RPKM to normalize the signal of 19 features for each exon/intron events from all samples
#RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 ) (from: https://www.metagenomics.wiki/pdf/qc/RPKM)
#geneLength is replaced by the exon or intron length

#import the count of reads of all data
epigenetic_feature_bam_read_count <- 
  read_csv("data/Chip_seq/no_treatment/bam/epigenetic_feature_bam_read_count.csv", col_names = c("epimark", "read_count"))

#join the read count of each epimark sample to tables
combined_all_refined_epimark <-
  bind_rows(refined_epimark) %>%
  left_join(epigenetic_feature_bam_read_count)

#calculate RPKM-like value of each signal for each event having peak signals 
#events with no peak signal will be included later
RI_HS0_PSI_epimark_RPKM_no_NA <-
  join_overlap_intersect(as_granges(RI_HS0_PSI), as_granges(combined_all_refined_epimark)) %>%
  as_tibble() %>%
  mutate(RPKM = fold_enrichment/(width/1000 * read_count/1000000)) %>%
  dplyr::select(seqnames, start, end, gene_name, AS_type, strand, PSI, PSI_group, epimark, RPKM)
  
#make the final raw table of all events with or without peak signals
#a trick to repeat dataframes 20 times
#https://www.geeksforgeeks.org/repeat-rows-of-dataframe-n-times-in-r/

RI_HS0_PSI_epimark_RPKM_all_events_raw <-
  map_dfr(seq_len(20), ~RI_HS0_PSI) %>%
  mutate(epimark = rep(epimark_HS0, each = nrow(RI_HS0_PSI))) %>%
  group_by(epimark, PSI_group) %>%
  mutate(event_label = c(1:n())) %>%
  mutate(PSI_group_event = str_c(PSI_group, "_", event_label)) %>%
  left_join(RI_HS0_PSI_epimark_RPKM_no_NA) %>%
  replace_na(list(RPKM = 0)) %>%
  ungroup()

write_delim(RI_HS0_PSI_epimark_RPKM_all_events_raw, "./data/refined_AS_bed_for_ML_HS0/RI_HS0_PSI_epimark_RPKM_all_events_raw.txt", delim = "\t", col_names = TRUE)

#To run boruta using events from PSI groups with similar expression levels
#I use the data of expression levels from HS0 for the normalization
#import groups of expression levels
RI_ctrl_AS_DAS_boxplot_bed6_output_raw <-
  read_delim("./data/M82_annotation_data/RI_ctrl_AS_DAS_gene_expression_level_0H_bed6_output.bed", delim = "\t", col_names = TRUE) %>%
  dplyr::select(gene_name, mean_0H, quantile_group_0H)

#check how many events will be kept if using 3 groups of expression levels
RI_HS0_PSI_epimark_RPKM_all_events_raw %>%
  filter(epimark == "Pol2") %>%
  left_join(RI_ctrl_AS_DAS_boxplot_bed6_output_raw) %>%
  drop_na() %>%
  mutate(expression_level_groups = if_else(mean_0H <= quantile(mean_0H, 0.33), "low",
                                           if_else(mean_0H <= quantile(mean_0H, 0.67) & mean_0H > quantile(mean_0H, 0.33), "medium", "high"))) %>%
  group_by(PSI_group, expression_level_groups) %>%
  summarise(count = n(), median_ex_l = median(mean_0H)) %>%
  View()

#make the raw table necessary for running boruta
RI_HS0_PSI_epimark_RPKM_all_events_ex_l_groups <-
  RI_HS0_PSI_epimark_RPKM_all_events_raw %>%
  left_join(RI_ctrl_AS_DAS_boxplot_bed6_output_raw) %>%
  drop_na() %>%
  mutate(expression_level_groups = if_else(mean_0H <= quantile(mean_0H, 0.33), "low",
                                           if_else(mean_0H <= quantile(mean_0H, 0.67) & mean_0H > quantile(mean_0H, 0.33), "medium", "high"))) %>%
  dplyr::select(PSI_group, expression_level_groups, PSI_group_event, epimark, RPKM) %>%
  mutate(epimark = str_c(epimark, "_RPKM"),
         PSI_group = factor(PSI_group, levels = str_c("PSI_", c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95")))) %>%
  spread(key = epimark, value = RPKM) %>%
  dplyr::select(-PSI_group_event) 


#Run boruta for pairwise comparisons of 7 PSI groups
HS0_expression_level_groups <- c("low", "medium", "high")

PSI_type <- str_c("PSI_", c("5", "5_20", "20_40", "40_60", "60_80", "80_95", "95"))

comb_pair_PSI_type <- combs(PSI_type, 2)

comb_pair_PSI_type

for (i in seq_along(HS0_expression_level_groups)) {
  for(j in seq_along(comb_pair_PSI_type[,1])) {
    table_for_boruta <-
      RI_HS0_PSI_epimark_RPKM_all_events_ex_l_groups %>%
      filter(expression_level_groups == HS0_expression_level_groups[[i]]) %>%
      filter(PSI_group== comb_pair_PSI_type[j, 1] | PSI_group == comb_pair_PSI_type[j, 2]) %>%
      dplyr::select(-expression_level_groups)
    
    RI_pairwise_PSI_ML_result <- Boruta(PSI_group~., table_for_boruta, doTrace = 2)
    
    data_RI <- RI_pairwise_PSI_ML_result$ImpHistory %>%
      as_tibble()
    
    write_delim(data_RI, str_c("./analysis/AS_HS0_RI_SE_ML_results/refined_PSI_comparisons/refined_RI_pairwise_comparison_", HS0_expression_level_groups[[i]], "_ex_l_", comb_pair_PSI_type[j,1], "_", comb_pair_PSI_type[j,2], "_boruta_Zscore"), col_names = TRUE, delim = "\t")
}
}

#tidy the boruta result and make plots
path_output_boruta_Zscore <-
  str_c("./analysis/AS_HS0_RI_SE_ML_results/refined_PSI_comparisons/refined_RI_pairwise_comparison_", rep(HS0_expression_level_groups, each = nrow(comb_pair_PSI_type)), "_ex_l_", rep(str_c(comb_pair_PSI_type[,1], "_", comb_pair_PSI_type[,2]), length(HS0_expression_level_groups)), "_boruta_Zscore")

pairwise_PSI_comp_Zscore_boruta <-
  map(path_output_boruta_Zscore, read_delim, delim = "\t", col_names = TRUE)

#create the function for having titles of ggplot
make_ggtitles <-
  function(a){
    ggtitle_raw_1 <-
      str_remove(a[1], "./analysis/AS_HS0_RI_SE_ML_results/refined_PSI_comparisons/refined_RI_pairwise_comparison_")
    
    ggtitle_raw_2 <-
      str_replace(ggtitle_raw_1, "_ex_l_", " expression levels at HS0 \n pairwise comparison: ")
    
    ggtitle_raw_3 <-
      str_remove(ggtitle_raw_2, "_boruta_Zscore")
    
    ggtitle_raw_4 <-
      str_replace(ggtitle_raw_3, "_PSI", " and PSI")
    
    str_c("Zscore of boruta for ", ggtitle_raw_4)
    
  }

#test
make_ggtitles(path_output_boruta_Zscore[[3]])

#create the function to make ggplot
tidy_boruta_table_and_make_plot <- function(data, a){
  ggtitle_output <- make_ggtitles(a[1])
  
  data %>%
    as_tibble() %>%
    pivot_longer(everything()) %>%
    filter(is.finite(value)==TRUE) %>%
    mutate(name = str_remove(name, "_RPKM")) %>%
    mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
    ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
    geom_boxplot(aes(fill = mark)) +
    ggtitle(ggtitle_output) +
    theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 36),
          legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
          axis.text.y = element_text(color = "black", size = 28, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
          axis.text.x = element_text(colour = "black", size = 32, face = "bold", angle = 45, vjust = 0.5))
}

ggplot_pairwise_PSI_comp_Zscore_boruta <-
  pmap(list(pairwise_PSI_comp_Zscore_boruta, path_output_boruta_Zscore), tidy_boruta_table_and_make_plot)

tidy_boruta_table_and_make_plot(pairwise_PSI_comp_Zscore_boruta[[1]], path_output_boruta_Zscore[[1]])

#create function to save ggplots
ggsave(str_c(path_output_boruta_Zscore[[1]], ".jpeg"), ggplot_pairwise_PSI_comp_Zscore_boruta[[1]], width = 400, height = 300, units = c("mm"), dpi = 320)

pwalk(list(str_c(path_output_boruta_Zscore, ".jpeg"), ggplot_pairwise_PSI_comp_Zscore_boruta), ggsave, width = 400, height = 300, units = c("mm"), dpi = 320)

str_c(path_output_boruta_Zscore[[1]], ".jpeg")

make_ggtitles(path_output_boruta_Zscore[[3]])


#create RI events from different PSI groups normalized by expression levels at HS0 for making deeptools plots
#splicing site and gene bodies of RI events will be both created for making plots

M82_all_gene_coordinates <-
  M82_rMATs_anno_all_gene %>%
  dplyr::select(gene_name, start_gene = str, end_gene = end)

RI_HS0_PSI_groups_normed_ex_l <-
  RI_HS0_PSI_epimark_RPKM_all_events_raw %>%
  filter(epimark == "Pol2") %>%
  left_join(RI_ctrl_AS_DAS_boxplot_bed6_output_raw) %>%
  drop_na() %>%
  mutate(expression_level_groups = if_else(mean_0H <= quantile(mean_0H, 0.33), "low",
                                           if_else(mean_0H <= quantile(mean_0H, 0.67) & mean_0H > quantile(mean_0H, 0.33), "medium", "high"))) %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, PSI_group, expression_level_groups) %>%
  left_join(M82_all_gene_coordinates)

create_output_normed_ex_l_PSI_groups <-
  function(a, b){
    output_PSI_group_ss <-
      RI_HS0_PSI_groups_normed_ex_l %>%
      filter(PSI_group == a[1] & expression_level_groups == b[1]) %>%
      dplyr::select(-start_gene, -end_gene)
    
    output_PSI_group_gene_body <-
    RI_HS0_PSI_groups_normed_ex_l %>%
      filter(PSI_group == a[1] & expression_level_groups == b[1]) %>%
      dplyr::select(seqnames, start = start_gene, end = end_gene, AS_type, gene_name, strand, PSI_group, expression_level_groups) %>%
      distinct()
    
    write_delim(output_PSI_group_gene_body, str_c("./data/refined_AS_bed_for_ML_HS0/refined_RI_gene_body_HS0_", a[1], "_", b[1], "_expression_levels_bed6.bed"), delim = "\t", col_names = FALSE)
    
    write_delim(output_PSI_group_ss, str_c("./data/refined_AS_bed_for_ML_HS0/refined_RI_ss_HS0_", a[1], "_", b[1], "_expression_levels_bed6.bed"), delim = "\t", col_names = FALSE)
  
}

walk2(rep(PSI_type, length(HS0_expression_level_groups)), rep(HS0_expression_level_groups, each = length(PSI_type)), create_output_normed_ex_l_PSI_groups)
    