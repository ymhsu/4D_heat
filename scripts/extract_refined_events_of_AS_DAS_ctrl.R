Packages <- c("plyranges", "tidyverse", "doParallel", "foreach", "nullranges", "caTools", "ggpubr", "fs", "HelloRanges")
lapply(Packages, library, character.only = TRUE)

#The main objective is to generate 3 alternative-splicing (AS) groups based on refined criteria.
#The refined 3 AS groups will be comparable for their profiles of epigenetic features.
#There are several refinements for 3 groups, mainly for ctrl and stable-PSI genes
#First, rMATs was used to separately detect AS from 3 heat conditions (24°C, 40°C for 1 h, and 40°C for 6h, they are below described as HS0/HS1/HS6)
#Second, the gene expression level of 3 AS groups will be normalized by the gene expression level of total expressed genes
#Gene sizes of 3 AS groups need to be concerned, and necessary refinement will be done 

#Here is the definition of 3 AS groups
#ctrl: genes without any AS (meaning they are constitutively spliced)
#stable PSI: genes which don't have differential spliced events between treatments
#DAS: gene with differential spliced events between treatments


##Define the group of ctrl and stable-PSI genes
#first, I focus on AS events found in 3 heat conditions (24°C, 40°C for 1 h, and 40°C for 6h)
#import AS events from 3 conditions (all of these event produced by rMATS 4.1.2)
#assign the variable for 4 types of AS
AS_type_import <- c("RI", "SE", "A5SS", "A3SS")

#create AS data with the average of PSI (since there are two replicates of RNA-seq data)

#produce the function to import data, tidy data, and select necessary columns for AS events from 3 heat conditions
produce_AS_3_trt <- function(a) {
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
  
  #set the path of AS imported from 3 directory (3 trt)
  path_input_AS_raw <-
    str_c("./data/AS_df_trt/HS", c(0, 1, 6), "/", a[1], ".MATS.JC.txt")
  
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
    dplyr::select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, GeneID) %>%
    mutate(GeneID = str_remove(GeneID, "gene:")) %>%
    mutate(AS_type = a[1]) %>%
    mutate(trt = c(
      rep("HS0", nrow_AS_df_trt[[1]]),
      rep("HS1", nrow_AS_df_trt[[2]]),
      rep("HS6", nrow_AS_df_trt[[3]])
    ))
}

#import data
all_AS_3_trt_separated_raw <-
  map(AS_type_import, produce_AS_3_trt) %>%
  bind_rows()

#calculate PSI of AS events 
all_AS_3_trt_separated_segments <- 
  all_AS_3_trt_separated_raw %>%
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
  dplyr::select(seqnames = chr, start = str, end, AS_type, gene_name = GeneID, strand, trt) %>%
  arrange(AS_type, trt, seqnames, start)

write_delim(all_AS_3_trt_separated_segments, "./data/AS_df_trt/all_AS_3_trt_combined_segments.bed", delim = "\t", col_names = TRUE)


#import the file of refined DAS genes
All_DAS_events_all_comparisons_refined_for_granges <-
  read_delim("./data/rMATS_out/All_DAS_events_all_comparisons_refined_for_granges", delim = "\t", col_names = TRUE)

#Produce the table with DAS and stable-PSI events
all_AS_DAS_combined_table <- 
  all_AS_3_trt_separated_segments %>%
  left_join(All_DAS_events_all_comparisons_refined_for_granges) %>%
  mutate(group_type = if_else(is.na(comp)==TRUE, "stable_PSI", "DAS"),
         comp = if_else(is.na(comp) == TRUE, "no", comp)) %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, trt, comp, group_type) 

#Here I wanted to check all DAS can be found in all AS separately detected by rMATS
#We should be able to find all events of DAS in the AS set 
All_DAS_events_all_comparisons_refined_for_granges %>%
  left_join(all_AS_3_trt_separated_segments) %>%
  filter(is.na(more_AS)==FALSE)

all_AS_DAS_combined_table %>%
  filter(gene_name == "Solyc02g079260.2")

#produce genes containing AS (stable PSI or DAS)
M82_all_genes_overlapped_AS <- 
  all_AS_DAS_combined_table %>%
  mutate(overlap_AS = "yes") %>%
  dplyr::select(gene_name, overlap_AS) %>%
  distinct()


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

#produce the raw table of introns with all DAS events (3 comp), AS events (events separately found in 3 trt), and ctrl genes larger than 3.5 kb(no AS)
#This refinement for gene size is because we found control genes having smaller size than AS-included genes
#And this difference leads to different profiles of epigenetic features
M82_all_gene_size <-
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name")) %>%
  mutate(gene_size = end - str) %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  filter(source == "Liftoff") %>%
  dplyr::select(seqnames = chr, start = str, end, type = feature, source, gene_name, strand) %>%
  mutate(size = end - start) %>%
  dplyr::select(gene_name, size)

#produce all introns of 3 AS groups from genes larger than 3.5 kb (including all types of AS)
all_AS_DAS_ctrl_introns_raw <-
  M82_rMATs_anno_all_intron %>%
  left_join(M82_all_genes_overlapped_AS) %>%
  replace_na(list(overlap_AS = "no")) %>%
  filter(overlap_AS != "yes") %>%
  mutate(AS_type = "noAS", trt = "all", comp = "no", group_type = "ctrl") %>%
  dplyr::select(seqnames, start, end, AS_type, gene_name, strand, trt, comp, group_type) %>%
  bind_rows(., all_AS_DAS_combined_table) %>%
  mutate(chr = str_remove(seqnames, "chr")) %>%
  filter(str_detect(chr, "[a-z]", negate = TRUE)) %>%
  mutate(chr = as.double(chr)) %>%
  drop_na() %>%
  arrange(group_type, AS_type, chr, start) %>%
  left_join(M82_all_gene_size) %>%
  filter(size > 3500) %>%
  dplyr::select(-size)


write_delim(all_AS_DAS_ctrl_introns_raw, "./data/AS_df_trt/all_AS_DAS_ctrl_introns_raw.bed", delim = "\t", col_names = TRUE)

#produce all introns of 3 AS groups from genes larger than 3.5 kb (including only RI, the comp between HS1 and HS6 is removed from DAS)
all_RI_stable_PSI_DAS_ctrl_introns_raw <-
  all_AS_DAS_ctrl_introns_raw %>%
  filter(AS_type == "RI" | AS_type == "noAS") %>%
  filter(comp != "HS1_HS6") %>%
  dplyr::select(-trt, -comp) %>%
  distinct() 

write_delim(all_RI_stable_PSI_DAS_ctrl_introns_raw, "./data/AS_df_trt/all_RI_stable_PSI_DAS_ctrl_introns_splicing_sites_raw.bed", delim = "\t", col_names = TRUE)

#second, the above-mentioned part gets focusing on each event representing splicing sites,
#we can create genes of DAS, stable PSI and ctrl focusing on gene bodies
#while we focused on gene bodies, the criterion for a gene being DAS is that any of all introns once recognized as a DAS
#Thus, we need to extract all genes ID having DAS introns, and remove introns as stable-PSI but with the same gene ID
#finally, stable-PSI genes are kept because none of their AS-included introns are defined as DAS events

#create the gene list of DAS
RI_DAS_gene_list <-
  all_RI_stable_PSI_DAS_ctrl_introns_raw %>%
  filter(group_type == "DAS") %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  mutate(DAS_or_not = "yes")


#remove the DAS gene from stable-PSI genes
RI_stable_PSI_ctrl_gene_list <-
  all_RI_stable_PSI_DAS_ctrl_introns_raw %>%
  filter(group_type != "DAS") %>%
  left_join(RI_DAS_gene_list) %>%
  replace_na(list(DAS_or_not = "no")) %>%
  filter(DAS_or_not == "no") %>%
  dplyr::select(gene_name, DAS_or_not, group_type) %>%
  distinct()

#combine DAS gene list and ctrl/stable gene list
RI_3_AS_group_gene_list <-
  RI_DAS_gene_list %>%
  mutate(group_type = "DAS") %>%
  bind_rows(., RI_stable_PSI_ctrl_gene_list) 

#based on annotation data (mRNA/gene) and the gene list of 3 AS groups
#we can produce mRNA of 3 AS groups

#import annotation data of all genes
M82_rMATs_anno_all_gene <- 
  read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_gene.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name")) %>%
  mutate(gene_size = end - str) %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  filter(source == "Liftoff") %>%
  dplyr::select(seqnames = chr, start = str, end, type = feature, source, gene_name, strand) 

#combine the gene list of 3 AS groups and the annotation data to get their regions 
RI_3_AS_group_mRNA <-
  M82_rMATs_anno_all_gene %>%
  left_join(RI_3_AS_group_gene_list) %>%
  drop_na() %>%
  dplyr::select(seqnames, start, end, strand, gene_name, group_type)

write_delim(RI_3_AS_group_mRNA, "./data/AS_df_trt/all_RI_stable_PSI_DAS_ctrl_introns_mRNA_raw.bed", delim = "\t", col_names = TRUE)

#use the file of RNA expression levels of all genes to check the expression level of three groups
#and classify 3 AS groups into subgroups depending on their expression levels
#the final result will be the comparison of 3 AS subgroups in which gene expression levels are similar 
#import mRNA of 3 groups (DAS, stable PSI and ctrl)
RI_3_AS_group_mRNA <-
  read_delim("./data/AS_df_trt/all_RI_stable_PSI_DAS_ctrl_introns_mRNA_raw.bed", delim = "\t", col_names = TRUE)

#import splicing sites of 3 groups
all_RI_stable_PSI_DAS_ctrl_introns_raw <-
  read_delim("./data/AS_df_trt/all_RI_stable_PSI_DAS_ctrl_introns_splicing_sites_raw.bed", delim = "\t", col_names = TRUE)

#import RNA expression level files to compare expression levels of no-AS, AS (stable PSI) and DAS groups
expression_level_raw <- read_csv("./data/expression_level.csv", col_names = TRUE) %>%
  mutate(ID = str_remove(ID, "gene:"))

expression_level <- expression_level_raw %>%
  mutate(mean_0H = (RNA_M82_HS_0h_rep1_S4 + RNA_M82_HS_0h_rep2_S5)/2,
         mean_1H = (RNA_M82_HS_1h_rep1_S6 + RNA_M82_HS_1h_rep2_S7)/2,
         mean_6H = (RNA_M82_HS_6h_rep1_S8 + RNA_M82_HS_6h_rep2_S9)/2) %>%
  dplyr::select(gene_name = ID, mean_0H, mean_1H, mean_6H)

#produce boxplots of expression levels at HS0 for three RI groups (all events)
#create the function to make boxplots
boxplot_PSI_expression_level_f <-
  function(data){
    data %>%
      ggplot(aes(x=factor(group_type), y=mean_0H, group = group_type)) + 
      geom_boxplot(fill='#A4A4A4', color="black")+
      #coord_cartesian(ylim = c(0, 15000)) + 
      scale_x_discrete(name = "PSI groups") +
      scale_y_continuous(name = "normalized read counts") +
      theme_classic() +
      theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 36),
            legend.title = element_blank(), axis.title.y = element_text(color = "black", size = 28, face = "bold"), 
            axis.title.x = element_text(colour = "black", size = 28, face = "bold"), 
            axis.text.y = element_text(color = "black", size = 22, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
            axis.text.x = element_text(colour = "black", size = 22, face = "bold")) 
  }

#create the boxplot using all genes of 3 AS groups
#the result shows that control genes have much lower expression levels than the other AS groups
RI_3_AS_group_gene_expression_level_plot <-
  RI_3_AS_group_mRNA %>%
  left_join(expression_level) %>%
  drop_na() %>%
  boxplot_PSI_expression_level_f() +
  coord_cartesian(ylim = c(0, 15000))

#In order to refine 3 AS groups to have their gene expression levels comparable for the peak profiles
#I used the expression level data from HS0 to generate subgroups of 3 AS groups
#based on RNA expression levels from HS0, I grouped all genes into 10 quantiles and use it as a pre-filtration
#to filter out 3 RI groups with similar RNA expression levels for comparing their profiles of epigenetic marks
all_gene_expression_level_0_removed <-
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
                                                                                             if_else(quantile(mean_0H, 0.8) <= mean_0H & mean_0H < quantile(mean_0H, 0.9), "q9", "q10"))))))))))


#combine the information of quantiles of expression levels at HS0 and RI from 3 groups together 
#for latter producing boxplots and bed6 files
#all genes having no expression at HS0 are removed (this could be changed while in other HS conditions)
RI_ctrl_AS_DAS_boxplot_bed6_output_raw <-
  RI_3_AS_group_mRNA %>%
  left_join(expression_level) %>%
  drop_na() %>%
  left_join(all_gene_expression_level_0_removed) %>%
  drop_na() %>%
  mutate(RI_group = if_else(group_type == "ctrl", "RI_ctrl_no_AS_genes", 
                            if_else(group_type == "stable_PSI", "RI_all_stable_PSI_genes", "RI_DAS_genes_H16_aft_H0"))) %>%
  dplyr::select(seqnames, start, end, gene_name, mean_0H, strand, group_type, RI_group, quantile_group_0H, mean_1H, mean_6H)


#check how many genes do we have for each quantile of 3 RI groups, and the statistics of expression levels for each quantile of 3 AS groups
RI_ctrl_AS_DAS_boxplot_bed6_output_raw %>%
  group_by(group_type, quantile_group_0H) %>%
  summarise(count = n(), median_ex_0H = median(mean_0H), median_ex_1H = median(mean_1H), median_ex_6H = median(mean_6H)) %>%
  View()

#make boxplots to check whether expression levels of 3 AS groups in each quantile are similar
RI_ctrl_AS_DAS_boxplot_bed6_output_raw %>%
  filter(quantile_group_0H != "q1" & quantile_group_0H != "q2" & quantile_group_0H != "q3") %>%
  boxplot_PSI_expression_level_f() +
  facet_wrap(~quantile_group_0H, scale = "free_y")

#produce files for quantiles by event numbers
produce_norm_expression_level_RI_AS_DAS <-
  function(a){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_raw %>%
      filter(quantile_group_0H != "q1" & quantile_group_0H != "q2") %>%
      filter(quantile_group_0H == a[1]) %>%
      split(.$RI_group)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level/", names(output_l), "_norm_expression_level_", output_l[[1]]$quantile_group_0H[[1]], "_ctrl_stable_PSI_corrected_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

walk(str_c("q", seq(3, 10)), produce_norm_expression_level_RI_AS_DAS)

#after modifying 3 AS groups, here we compare the content of them
#we didn't check the content in each quantile, we combined all quantiles of each AS groups for the comparison
#the idea is identify whether it is comparable for these expressed genes (the first criterion) with gene size larger than 3.5 (the second criterion) from 3 AS groups 
#we then check gene sizes, intron numbers, percentage occupied by introns, and intron sizes between 3 AS groups

#calculate the simple statistics of introns belonging to genes in annotation data
M82_rMATs_anno_all_intron_stat <-
  M82_rMATs_anno_all_intron_exon %>%
  group_by(gene_name, feature) %>%
  summarise(count_feature = n(), sum_bp_feature = sum(end - start)) %>%
  group_by(gene_name) %>%
  mutate(sum_bp_total = sum(sum_bp_feature)) %>%
  mutate(proportion_feature = sum_bp_feature/sum_bp_total) %>%
  filter(feature == "intron") %>%
  ungroup()

#create function for modifying ggplot
theme_ym_intron_stat <-
  function(data) {
    data +
      ylab("proportion in its subgroup") +
      theme_classic() +
      theme(strip.text.x = element_text(colour = "black", face = "bold", size = 28), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 36),
            legend.title = element_blank(), axis.title.y = element_text(color = "black", size = 28, face = "bold"), 
            axis.title.x = element_text(colour = "black", size = 28, face = "bold"), 
            axis.text.y = element_text(color = "black", size = 22, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
            axis.text.x = element_text(colour = "black", size = 18, face = "bold"))
  }

#gene size
rawplot_gene_size_3_AS_group <-
  RI_3_AS_group_mRNA %>%
  left_join(all_gene_expression_level_0_removed) %>%
  drop_na() %>%
  ggplot(aes(x = (end - start)/1000)) +
  geom_histogram(aes(y = ..density..), binwidth=1) +
  facet_wrap(~group_type, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  coord_cartesian(xlim=c(0, 40)) +
  xlab("The size of genes (kb)")

fplot_gene_size_3_AS_group <-
  theme_ym_intron_stat(rawplot_gene_size_3_AS_group)

ggsave("./analysis/The_statistics_of_introns_among_3_AS_group/gene_size_3_AS_group.jpeg", fplot_gene_size_3_AS_group, width = 400, height = 240, units = c("mm"), dpi = 320)

#count of introns
rawplot_count_intron_RI_3_AS_group_gene <- 
  RI_3_AS_group_mRNA %>%
  left_join(all_gene_expression_level_0_removed) %>%
  left_join(M82_rMATs_anno_all_intron_stat) %>%
  drop_na() %>%
  ggplot(aes(x = count_feature)) +
  geom_histogram(aes(y = ..density..), binwidth=1) +
  facet_wrap(~group_type) +
  scale_x_continuous(breaks = seq(0, 40, 5)) +
  coord_cartesian(xlim=c(0, 40)) +
  xlab("Number of introns in genes") 

fplot_count_intron_RI_3_AS_group_gene <-
  theme_ym_intron_stat(rawplot_count_intron_RI_3_AS_group_gene)

ggsave("./analysis/The_statistics_of_introns_among_3_AS_group/count_intron_RI_3_AS_group_gene.jpeg", fplot_count_intron_RI_3_AS_group_gene, width = 400, height = 240, units = c("mm"), dpi = 320)

#proportion of introns in each gene
rawplot_proportion_intron_RI_3_AS_group_gene <- 
  RI_3_AS_group_mRNA %>%
  left_join(all_gene_expression_level_0_removed) %>%
  left_join(M82_rMATs_anno_all_intron_stat) %>%
  drop_na() %>%
  ggplot(aes(x = proportion_feature)) +
  geom_histogram(aes(y = 0.1*..density..), binwidth=0.1) +
  facet_wrap(~group_type) +
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  xlab("The proportion of introns in genes")

fplot_proportion_intron_RI_3_AS_group_gene <-
  theme_ym_intron_stat(rawplot_proportion_intron_RI_3_AS_group_gene)


ggsave("./analysis/The_statistics_of_introns_among_3_AS_group/proportion_intron_RI_3_AS_group_gene.jpeg", fplot_proportion_intron_RI_3_AS_group_gene, width = 400, height = 240, units = c("mm"), dpi = 320)

#intron size
list_3_AS_group_mRNA <-
  RI_3_AS_group_mRNA %>%
  left_join(all_gene_expression_level_0_removed) %>%
  drop_na() %>%
  dplyr::select(gene_name, group_type)

rawplot_intron_size_3_AS_group_gene <-
  M82_rMATs_anno_all_intron_exon %>%
  filter(feature == "intron") %>%
  left_join(list_3_AS_group_mRNA) %>%
  drop_na() %>%
  mutate(intron_size = end - start) %>%
  ggplot(aes(x = intron_size)) +
  geom_histogram(aes(y = 100*..density..),
                 breaks = seq(from = 0, to = 1000, by = 100), binwidth=100) +
  facet_wrap(~group_type, scales = "free_x") +
  scale_x_continuous(breaks = seq(0, 1000, 200)) +
  coord_cartesian(xlim=c(0, 1000)) +
  xlab("The size of introns (bp)")

fplot_intron_size_3_AS_group_gene <-
  theme_ym_intron_stat(rawplot_intron_size_3_AS_group_gene)

ggsave("./analysis/The_statistics_of_introns_among_3_AS_group/intron_size_3_AS_group_gene.jpeg", fplot_intron_size_3_AS_group_gene, width = 400, height = 240, units = c("mm"), dpi = 320)

#after creating gene bodies of all AS events from 3 AS groups
#we created the splicing site of each event

#create the list of gene with their quantiles based on HS0 gene expression levels
all_gene_expression_level_0_removed_light <-
  all_gene_expression_level_0_removed %>%
  dplyr::select(gene_name, quantile_group_0H)


RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site <-
  all_RI_stable_PSI_DAS_ctrl_introns_raw %>%
  left_join(expression_level) %>%
  drop_na() %>%
  left_join(all_gene_expression_level_0_removed_light) %>%
  drop_na() %>%
  mutate(RI_group = if_else(group_type == "ctrl", "RI_ctrl_no_AS_genes", 
                            if_else(group_type == "stable_PSI", "RI_all_stable_PSI_genes", "RI_DAS_genes_H16_aft_H0"))) %>%
  dplyr::select(seqnames, start, end, gene_name, mean_0H, strand, group_type, RI_group, quantile_group_0H, mean_1H, mean_6H)


RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site %>%
  group_by(group_type, quantile_group_0H) %>%
  summarise(count = n()) %>%
  View()

RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site %>%
  #filter(quantile_group_0H != "q1" & quantile_group_0H != "q2" & quantile_group_0H != "q3" & quantile_group_0H != "q4"& quantile_group_0H != "q5") %>%
  group_by(group_type, quantile_group_0H) %>%
  summarise(count = n(), median_ex = median(mean_0H)) %>%
  View()
boxplot_PSI_expression_level_f() +
  #ylim(200, 6000) +
  facet_wrap(~quantile_group_0H, scale = "free_y")

#produce these events in 10 groups
#ss means splicing sites
produce_norm_expression_level_RI_AS_DAS_ss <-
  function(a){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site %>%
      filter(quantile_group_0H != "q1" & quantile_group_0H != "q2") %>%
      filter(quantile_group_0H == a[1]) %>%
      split(.$RI_group)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level/", names(output_l), "_norm_expression_level_", output_l[[1]]$quantile_group_0H[[1]], "_ctrl_stable_PSI_corrected_ss_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

walk(str_c("q", seq(3, 10)), produce_norm_expression_level_RI_AS_DAS_ss)


RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site %>%
  filter(quantile_group_0H != "q1" & quantile_group_0H != "q2") %>%
  group_by(group_type, quantile_group_0H) %>%
  summarise(median_size = median(end-start)) %>%
  View()



#Jeremie proposed another way to establish groups in terms of expression levels
#which is setting the range of expression levels
#check the boundaries of the minimum and maximum expression levels
#produce the gene bodies including 3 AS-group events
expression_level %>%
  filter(mean_0H > 5 & mean_0H < 50000) %>%
  ggplot() +
  geom_histogram(aes(log(mean_0H)), bins = 150)

boundary_expression_level <- seq(log(5), log(50000), length.out = 11)

all_gene_expression_level_log_group_0_removed <-
  M82_rMATs_anno_all_gene %>%
  left_join(expression_level) %>%
  drop_na() %>%
  filter(mean_0H > 5 & mean_0H < 50000) %>%
  mutate(log_mean_0H = log(mean_0H)) %>%
  #View()
  mutate(quantile_group_0H = if_else(log_mean_0H > boundary_expression_level[[1]] & log_mean_0H <= boundary_expression_level[[2]], "q1", 
                                   if_else(log_mean_0H > boundary_expression_level[[2]] & log_mean_0H <= boundary_expression_level[[3]], "q2",
                                           if_else(log_mean_0H > boundary_expression_level[[3]] & log_mean_0H <= boundary_expression_level[[4]], "q3",
                                                   if_else(log_mean_0H > boundary_expression_level[[4]] & log_mean_0H <= boundary_expression_level[[5]], "q4",
                                                           if_else(log_mean_0H > boundary_expression_level[[5]] & log_mean_0H <= boundary_expression_level[[6]], "q5",
                                                                   if_else(log_mean_0H > boundary_expression_level[[6]] & log_mean_0H <= boundary_expression_level[[7]], "q6",
                                                                           if_else(log_mean_0H > boundary_expression_level[[7]] & log_mean_0H <= boundary_expression_level[[8]], "q7",
                                                                                   if_else(log_mean_0H > boundary_expression_level[[8]] & log_mean_0H <= boundary_expression_level[[9]], "q8",
                                                                                           if_else(log_mean_0H > boundary_expression_level[[9]] & log_mean_0H <= boundary_expression_level[[10]], "q9", "q10")))))))))) 
                                            

RI_ctrl_AS_DAS_boxplot_bed6_output_raw_ex_log <-
  RI_3_AS_group_mRNA %>%
  left_join(expression_level) %>%
  drop_na() %>%
  left_join(all_gene_expression_level_log_group_0_removed) %>%
  drop_na() %>%
  mutate(RI_group = if_else(group_type == "ctrl", "RI_ctrl_no_AS_genes", 
                            if_else(group_type == "stable_PSI", "RI_all_stable_PSI_genes", "RI_DAS_genes_H16_aft_H0"))) %>%
  dplyr::select(seqnames, start, end, gene_name, mean_0H, strand, group_type, RI_group, quantile_group_0H, mean_1H, mean_6H)


RI_ctrl_AS_DAS_boxplot_bed6_output_raw_ex_log %>%
  group_by(RI_group, quantile_group_0H) %>%
  summarise(count = n()) %>%
  View()


#produce files for quantiles by the range of expression levels
produce_norm_expression_level_range_RI_AS_DAS <-
  function(a){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_raw_ex_log %>%
      filter(quantile_group_0H == a[1]) %>%
      split(.$RI_group)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level/", names(output_l), "_norm_expression_level_range_", output_l[[1]]$quantile_group_0H[[1]], "_ctrl_stable_PSI_corrected_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

map(str_c("q", seq(1, 10)), produce_norm_expression_level_range_RI_AS_DAS)

#produce splicing sites of 3 AS-group events
all_gene_expression_level_log_group_0_removed_light <-
  all_gene_expression_level_log_group_0_removed %>%
  dplyr::select(gene_name, quantile_group_0H)

#produce files based on quantiles of ranges of expression levels
all_gene_expression_level_log_group_0_removed_light <-
  all_gene_expression_level_log_group_0_removed %>%
  dplyr::select(gene_name, quantile_group_0H)

RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site_ex_l_range <-
  all_RI_stable_PSI_DAS_ctrl_introns_raw %>%
  left_join(expression_level) %>%
  drop_na() %>%
  left_join(all_gene_expression_level_log_group_0_removed_light) %>%
  drop_na() %>%
  mutate(RI_group = if_else(group_type == "ctrl", "RI_ctrl_no_AS_genes", 
                            if_else(group_type == "stable_PSI", "RI_all_stable_PSI_genes", "RI_DAS_genes_H16_aft_H0"))) %>%
  dplyr::select(seqnames, start, end, gene_name, mean_0H, strand, group_type, RI_group, quantile_group_0H, mean_1H, mean_6H)


RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site_ex_l_range %>%
  group_by(group_type, quantile_group_0H) %>%
  summarise(count = n()) %>%
  View()

RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site_ex_l_range %>%
  #filter(quantile_group_0H != "q1" & quantile_group_0H != "q2" & quantile_group_0H != "q3" & quantile_group_0H != "q4"& quantile_group_0H != "q5") %>%
  #group_by(group_type, quantile_group_0H) %>%
  #summarise(count = n(), median_ex = median(mean_0H)) %>%
  #View()
  boxplot_PSI_expression_level_f() +
  #ylim(200, 6000) +
  facet_wrap(~quantile_group_0H, scale = "free_y")


#ss means splicing sites
produce_norm_expression_level_range_RI_AS_DAS_ss <-
  function(a){
    
    output_l <-
      RI_ctrl_AS_DAS_boxplot_bed6_output_raw_splicing_site_ex_l_range %>%
      filter(quantile_group_0H == a[1]) %>%
      split(.$RI_group)
    
    #output_l[[1]]$quantile_group_0H[[1]]
    names(output_l)
    
    output_name <- str_c("./data/ctrl_AS_DAS_genes_normed_expression_level/", names(output_l), "_norm_expression_level_range_", output_l[[1]]$quantile_group_0H[[1]], "_ctrl_stable_PSI_corrected_ss_bed6.bed")
    
    pmap(list(output_l, output_name), write_delim, delim ="\t", col_names = FALSE)
    
  }

walk(str_c("q", seq(1, 10)), produce_norm_expression_level_range_RI_AS_DAS_ss)











#make all Liftoff genes as gr file
M82_rMATs_anno_all_gene_gr <- 
  M82_rMATs_anno_all_gene %>%
  dplyr::select(seqnames = chr, start = str, end, strand, feature, source, gene_name) %>%
  filter(source == "Liftoff") %>%
  mutate(gene_name = str_replace(gene_name, "(.*\\..*)\\..*", "\\1")) %>%
  as_granges

all_AS_HS0_gene_name <- 
  join_overlap_intersect(M82_rMATs_anno_all_gene_gr, as_granges(all_AS_HS0_f)) %>%
  as_tibble() %>%
  dplyr::select(gene_name, AS_type) %>%
  distinct()

test <-
  join_overlap_intersect(M82_rMATs_anno_all_gene_gr, as_granges(All_DAS_events_all_comparisons_refined_for_granges)) %>%
  as_tibble() %>%
  distinct()

all_AS_HS0_gene_name %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  left_join(test)

#test HS0 
col_names_uni <- c("ID_1", "GeneID", "geneSymbol", "chr", "strand", "pos_1", "pos_2", "pos_3", "pos_4", "pos_5", 
                   "pos_6", "ID_2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen",
                   "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")

RI_rMATS_4_1_2 <-
read_delim("./data/AS_df_trt/HS0/RI.MATS.JC.txt", col_names = col_names_uni, delim = "\t", skip = 1) %>%
  mutate(AS_type = "RI") %>%
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
  dplyr::select(seqnames = chr, start = str, end, AS_type) %>%
  mutate(rMATS_version_new = "new")


read_delim("./data/AS_df_trt/HS0/RI.MATS.JC.txt", col_names = col_names_uni, delim = "\t", skip = 1) %>%
  mutate(AS_type = "RI") %>%
  dplyr::select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, IncLevel1, AS_type) %>%
  filter(str_detect(IncLevel1, "NA")==TRUE)

RI_rMATS_4_1_0 <- 
read_delim("./data/AS_events_HS0/RI.MATS.JC.txt", col_names = col_names_uni, delim = "\t", skip = 1) %>%
  mutate(AS_type = "RI") %>%
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
  dplyr::select(seqnames = chr, start = str, end, AS_type) %>%
  mutate(rMATS_version_old = "old") 

RI_rMATS_4_1_2 %>%
  left_join(RI_rMATS_4_1_0) %>%
  filter(is.na(rMATS_version_old)==TRUE)

RI_rMATS_4_1_0 %>%
  left_join(RI_rMATS_4_1_2) %>%
  filter(is.na(rMATS_version_new)==TRUE)

16370 662

15748/16370
all_AS_HS0_f %>%
  filter(AS_type == "RI") %>%
  left_join(RI_rMATS_4_1_2) %>%
  filter(is.na(rMATS_version)==TRUE)
15742

#probably deleted#
#here, AS could happen in exons of a mRNA, so we need to use all mRNA as the anchor 
#to make the intersection between genes and AS/DAS, then get the gene set having always spliced introns
M82_Liftoff_all_mRNA <- 
  read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12_Chr_corrected_mRNA.bed", delim = "\t", 
             col_names = c("seqnames", "start", "end", "c4", "c5", "strand", "source", "type", "c7", "anno")) %>%
  mutate(anno = str_remove(anno, "ID=mRNA:")) %>%
  mutate(c7 = str_remove(anno, ";Name.*")) %>%
  mutate(c7 = str_replace(c7, "(.*\\..*)\\..*", "\\1")) %>%
  dplyr::select(seqnames, start, end, type, source, gene_name = c7, strand) 

test_mRNA <- 
  M82_Liftoff_all_mRNA %>%
  dplyr::select(chr = seqnames, str = start, end) %>%
  mutate(mRNA = "yes")

M82_rMATs_anno_all_gene %>%
  filter(source != "REGARN") %>%
  left_join(test_mRNA) %>%
  filter(is.na(mRNA)==TRUE)

gene_having_intron <-
  M82_rMATs_anno_all_intron %>%
  dplyr::select(gene_name) %>%
  distinct() %>%
  mutate(gene_with_intron = "yes")

M82_Liftoff_gene_with_intron_mRNA <-
  M82_Liftoff_all_mRNA %>%
  left_join(gene_having_intron) %>%
  drop_na() %>%
  dplyr::select(-gene_with_intron)

M82_Liftoff_gene_with_intron_mRNA %>%
  filter(gene_name == "Solyc02g079260.2")

Solyc02g079260.2
#probably deleted


#deleted
M82_all_genes_overlapped_AS <-
  join_overlap_intersect(as_granges(M82_Liftoff_gene_with_intron_mRNA), as_granges(all_AS_DAS_combined_table)) %>%
  as_tibble() %>%
  dplyr::select(seqnames, start, end, gene_name.x, strand) %>%
  distinct() %>%
  mutate(overlap_AS = "yes") %>%
  dplyr::select(gene_name = gene_name.x, overlap_AS) %>%
  distinct()
#deleted