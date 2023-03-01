#Focusing on AS events identified in rMATs only for HS0
#The following process is to produce 7 classes of PSI for all AS and their internal controls
Sys.setenv(LANG = "en")
library(tidyverse)

#The modification of the header for all output files of rMATs
col_names_uni <- c("ID_1", "GeneID", "geneSymbol", "chr", "strand", "pos_1", "pos_2", "pos_3", "pos_4", "pos_5", 
                   "pos_6", "ID_2", "IJC_SAMPLE_1", "SJC_SAMPLE_1", "IJC_SAMPLE_2", "SJC_SAMPLE_2", "IncFormLen",
                   "SkipFormLen", "PValue", "FDR", "IncLevel1", "IncLevel2", "IncLevelDifference")

#create the function of producing 
produce_AS_HS0_PSI <- function(a){
  #import raw data of all .MATS.JC.txt files
  #variable "a" is different types of AS
  data_raw <- read_delim(str_c("./data/AS_events_HS0/", a[1], ".MATS.JC.txt"), col_names = col_names_uni, skip = 1) %>%
    select(chr, strand, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6, IncLevel1) %>%
    filter(str_detect(IncLevel1, "NA")==FALSE) %>%
    mutate(AS_type = a[1])
  
  #split IncLevel1 (there are two replicates which will be used for calculating average PSI)
  IncLevel1_m <- str_split(data_raw$IncLevel1, ",", simplify = TRUE)
  
  #add PSI data from two replicates to compute "PSI" which is the average PSI of r1 and r2 
  data_raw %>%
    mutate(
      IncLevel_r1 = as.double(IncLevel1_m[,1]),
      IncLevel_r2 = as.double(IncLevel1_m[,2])
    ) %>%
    mutate(PSI = (IncLevel_r1 + IncLevel_r2)/2) %>%
    select(-IncLevel_r1, -IncLevel_r2)
}

#assign the variable for 4 types of AS
AS_type <- c("RI", "SE", "A5SS", "A3SS")

#create AS data with the average of PSI
AS_PSI_raw <- 
AS_type %>%
  map(., produce_AS_HS0_PSI) %>%
  bind_rows() %>%
  mutate(AS_type = if_else(AS_type == "A5SS" | AS_type == "A3SS", str_remove(AS_type, "S"), AS_type))


AS_PSI_raw %>%
  View()

AS_PSI_raw %>%
  group_by(AS_type) %>%
  summarise()

#produce the list of 4 types of AS with 7 classes of PSI (nested list with two layers)
#In brief, PSI_5: PSI <= 0.05, PSI_5_20:  0.05 < PSI <= 0.2, etc.
AS_HSO_PSI_l <- 
AS_PSI_raw %>%
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
  mutate(PI = if_else(PSI <= 0.05, "PSI_5",
                      if_else(PSI > 0.05 & PSI <= 0.2, "PSI_5_20",
                      if_else(PSI > 0.2 & PSI <= 0.4, "PSI_20_40", 
                              if_else(PSI > 0.4 & PSI <= 0.6, "PSI_40_60",
                                      if_else(PSI > 0.6 & PSI <= 0.8, "PSI_60_80", 
                                              if_else(PSI > 0.8 & PSI <= 0.95, "PSI_80_95", "PSI_95"))))))) %>%
  mutate(trt = "HS0") %>%
  select(chr, str, end, trt, AS_type, PI) %>%
  distinct() %>%
  split(.$AS_type) %>%
  map(. %>% split(.$PI))


#make the output files of AS with PSI for the following analysis
#For example, we use these regions to produce deeptool plots
for (i in seq_along(AS_HSO_PSI_l)) {
  for (j in seq_along(AS_HSO_PSI_l[[1]])) {
    output_path <- str_c("./data/AS_events_HS0/AS_HS0_PSI_segment/AS_HS0_", AS_HSO_PSI_l[[i]][[j]]$AS_type[[1]],
                         "_", AS_HSO_PSI_l[[i]][[j]]$PI[[1]], ".bed")
    
    write_delim(AS_HSO_PSI_l[[i]][[j]], output_path, col_names = FALSE, delim = "\t")
  }
}

bind_rows(AS_HSO_PSI_l[[1]]) %>%
  distinct()

#produce the list of AS only with one layer (only separate 4 types of AS into 4 files)
#meaning 7 classes of PSI were combined into one table of each AS type
#These files are used for the intersection with genomic features,
#and then define segments for machine learning
AS_HSO_PSI_l_one_layer <-
AS_HSO_PSI_l %>%
  map(. %>% bind_rows()) %>%
  map(. %>% arrange(chr, str))


for (i in seq_along(AS_HSO_PSI_l_one_layer)) {
  output_path <- str_c("./data/AS_events_HS0/AS_HS0_PSI_segment/AS_HS0_", AS_HSO_PSI_l_one_layer[[i]]$AS_type[[1]],".bed")
  
  write_delim(AS_HSO_PSI_l_one_layer[[i]], output_path, col_names = FALSE, delim = "\t")
}

#check the size of AS
AS_PSI_raw %>%
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
  mutate(PI = if_else(PSI <= 0.05, "PSI_5",
                      if_else(PSI > 0.05 & PSI <= 0.2, "PSI_5_20",
                              if_else(PSI > 0.2 & PSI <= 0.4, "PSI_20_40", 
                                      if_else(PSI > 0.4 & PSI <= 0.6, "PSI_40_60",
                                              if_else(PSI > 0.6 & PSI <= 0.8, "PSI_60_80", 
                                                      if_else(PSI > 0.8 & PSI <= 0.95, "PSI_80_95", "PSI_95"))))))) %>%
  select(chr, str, end, AS_type, PI) %>%
  arrange(chr, str) %>%
  group_by(AS_type, PI) %>%
  #filter(end-str<10)
  summarise(median_region = median(end-str), mean_region = mean(end-str), count = n()) %>%
  View()


#import all features (exon/intron) intersected with AS
AS_type = c("A5S", "A3S", "SE", "RI")
feature = c("exon", "intron")

454436
all_AS_intersected_feature_raw <- tibble()

for (i in seq_along(AS_type)) {
  for (j in seq_along(feature)) {
    path_AS_mark <- str_c("./data/Intersected_AS_HS0_M82_anno_rMATs/", "AS_HS0_", AS_type[[i]], "_intersected_M82_rMATs_", feature[[j]], "_order.bed")
    
    a <- read_delim(path_AS_mark, delim = "\t", 
                    col_names = c("chr", "str", "end", "trt", "AS", "PI", "chr_2", "str_2", "end_2", "strand", "feature", "source", "anno", "order")) %>%
      select(feature, anno, order, trt, AS, PI) %>%
      distinct()
    
    all_AS_intersected_feature_raw <- bind_rows(all_AS_intersected_feature_raw, a)
      }
    }

#combine all genes intersected with AS (four types)  
gene_having_AS <- all_AS_intersected_feature_raw %>%
  select(feature, anno) %>%
  arrange(anno) %>%
  mutate(having_AS = "yes") %>%
  distinct()

##Import the information of all genes and classify them based on PSI
##the following analysis was to produce genes with RI or SE from 7 PSI groups
##And check the global signals of different epigenetic features
M82_mRNA <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12_Chr_corrected_mRNA.bed",
           col_names = c("chr", "str", "end", "name", "score", "strand", "source", "feature", "col9", "info"))

#add the information of gene name into additional column
M82_mRNA_gene_name <- M82_mRNA %>%
  mutate(anno = str_remove(info, ";Name=.*")) %>%
  mutate(anno = str_remove(anno, "ID=mRNA:"))

#produce list of gene names of different AS type (RI/SE)
AS_PSI_gene_name <- all_AS_intersected_feature_raw %>%
  select(anno, PI, AS) %>%
  distinct() 

#create the function for having the output of separate bed files of PSIs of each AS type (SE/RI)
producing_output_all_mRNA_AS_PSI_f <- function(a){
  no_AS_PSI_table <- M82_mRNA_gene_name %>%
    left_join(AS_PSI_gene_name) %>%
    filter(is.na(PI)==TRUE) %>%
    replace_na(list(PI = "no_PSI", AS = a[1])) 
  
  AS_PSI_table <- M82_mRNA_gene_name %>%
    left_join(AS_PSI_gene_name) %>%
    filter(AS == a[1])
  
  data_l <- bind_rows(no_AS_PSI_table, AS_PSI_table) %>%
    arrange(chr, str) %>%
    ungroup() %>%
    split(.$PI)
  pwalk(list(data_l, str_c("./analysis/AS_RI_SE_PSI_all_genes_bed/M82_rMATs_all_mRNA_", a[1], "_", names(data_l), ".bed")),
      write_delim, col_names = FALSE, delim = "\t" )
}

producing_output_all_mRNA_AS_PSI_f("RI")
producing_output_all_mRNA_AS_PSI_f("SE")


#Import all exon and intron with the order into
M82_rMATs_anno_all_exon_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))
M82_rMATs_anno_all_intron_order <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron_order.bed", col_names = c("chr", "str", "end", "strand", "feature", "source", "anno", "order"))

##check the distribution of SE and RI occurrence along the gene bodies
#extract all AS events (SE/RI) from the AS whole dataset
all_AS_intersected_feature_raw_SE <- all_AS_intersected_feature_raw %>%
  filter(AS == "SE")

all_AS_intersected_feature_raw_RI <- all_AS_intersected_feature_raw %>%
  filter(AS == "RI")

#produce the list of all genes for the following join with dataset having AS events
M82_rMATs_anno_all_feature_min_max <- 
  bind_rows(M82_rMATs_anno_all_exon_order, M82_rMATs_anno_all_intron_order) %>%
  arrange(chr, str) %>%
  group_by(feature, anno) %>%
  mutate(min_feature_order = min(order), max_feature_order = max(order)) %>%
  ungroup() %>%
  select(feature, anno, min_feature_order, max_feature_order) %>%
  distinct() %>%
  split(.$feature)

#calculate relative position of features with AS and make plot 
#genes with at least 3 exon/introns
#produce list of genes with the information of the order of the first, last and AS-event contained feature
#SE
data_SE_combined <- M82_rMATs_anno_all_feature_min_max$exon %>%
  filter(max_feature_order >= 3) %>%
  left_join(all_AS_intersected_feature_raw_SE) %>%
  drop_na()
#RI
data_RI_combined <- M82_rMATs_anno_all_feature_min_max$intron %>%
  filter(max_feature_order >= 3) %>%
  left_join(all_AS_intersected_feature_raw_RI) %>%
  drop_na()

#combine the data of SE and RI and calculate the relative position of AS
#relative position os AS is to use the order of AS-event contained feature subtracting the order of the first feature
#divided by the order of the last feature subtracting the order of the first feature
output_plot_raw <- bind_rows(data_SE_combined, data_RI_combined) %>%
  mutate(relative_position_AS = (order-min_feature_order)/(max_feature_order-min_feature_order)) %>%
  ggplot() +
  geom_histogram(aes(relative_position_AS), bins = 51) +
  facet_wrap(~feature)

output_plot <- theme_ym(output_plot_raw)

ggsave(str_c("./analysis/AS_RI_SE_location_plots/", "AS_feature_relative_position_at_least_3_events", ".jpeg"), output_plot, width = 400, height = 240, units = c("mm"), dpi = 320)

#check the bias of the SE/AS occurrence along the gene bodies
#it turned out SE prefers to locate in 5 prime regions (2% more) and RI prefers to located in 3 prime regions (5% more)
bind_rows(data_SE_combined, data_RI_combined) %>%
  mutate(relative_position_AS = (order-min_feature_order)/(max_feature_order-min_feature_order)) %>%
  #remove the middle (0.5)
  filter(relative_position_AS != 0.5) %>%
  mutate(posi_two_part = if_else(relative_position_AS <= 0.5, "5_prime", "3_prime")) %>%
  group_by(feature, posi_two_part) %>%
  summarise(count = n())

#Separated plots for different numbers of features are also made to see the overall trend
for (i in seq(2,30)) {
  data_SE <- M82_rMATs_anno_all_exon_min_max %>%
    filter(max_feature_order == i) %>%
    left_join(all_AS_intersected_feature_raw_SE) %>%
    drop_na()
  
  data_RI <- M82_rMATs_anno_all_intron_min_max %>%
    filter(max_feature_order == i) %>%
    left_join(all_AS_intersected_feature_raw_RI) %>%
    drop_na()
  
  output_plot_raw_separated <- bind_rows(data_SE, data_RI) %>%
    mutate(relative_position_AS = (order-min_feature_order)/(max_feature_order-min_feature_order)) %>%
    ggplot() +
    geom_histogram(aes(relative_position_AS), bins = 51) +
    facet_wrap(~feature)
  
  output_plot_separated <- theme_ym(output_plot_raw_separated)
  
  
  ggsave(str_c("./analysis/AS_RI_SE_location_plots/AS_feature_relative_position_", i, ".jpeg"), output_plot_separated, width = 400, height = 240, units = c("mm"), dpi = 320)
  
}


#combine all exon and intron as a list
all_features_M82_rMATs_anno_l <- list(M82_rMATs_anno_all_exon_order, M82_rMATs_anno_all_intron_order)
names(all_features_M82_rMATs_anno_l) <- c("exon", "intron")


#gene with AS
gene_having_AS_info_l_f <- all_features_M82_rMATs_anno_l %>%
  map(. %>% arrange(chr, str)) %>%
  map(. %>% left_join(gene_having_AS)) %>%
  map(. %>% drop_na()) %>%
  map(. %>% select(-having_AS)) %>%
  map(. %>% distinct())


#create the function for deciding which gene should be kept for the following ML analysis
#1. remove promoter effects (the first/second exon/intron are removed)
list_for_ML_data_f <- function(data, a) {
  data %>%
    #remove genes with only two exons/introns
    group_by(anno) %>%
    mutate(feature_count = n()) %>%
    filter(feature_count > 2) %>%
    select(-feature_count) %>%
    ungroup() %>%
    left_join(all_AS_intersected_feature_raw, by = c("feature", "anno", "order")) %>%
    replace_na(list(trt = "no", AS = "no", PI = "no", event = "no")) %>%
    #remove events in the first or second exon/intron (promoter effects)
    filter(order != 1 & order != 2) %>%
    #extract the specific type of AS (SE or RI)
    filter(AS != "A5S" & AS != "A3S" & AS != a[1]) %>%
    mutate(AS_type_label = if_else(AS == a[2], 1, 0)) %>%
    group_by(anno) %>%
    mutate(sum_AS_type_label = sum(AS_type_label)) %>%
    filter(sum_AS_type_label != 0) %>%
    #filter out genes without having AS
    mutate(no_event = if_else(trt == "no" & AS == "no" & PI == "no", 0, 1)) %>%
    group_by(anno) %>%
    mutate(sum_no_event = sum(no_event)) %>%
    filter(sum_no_event != 0) %>%
    ungroup() %>%
    select(-no_event, -sum_no_event, -AS_type_label, -sum_AS_type_label) %>%
    split(.$anno)
}

#create the list genes with sig events for the random selection of event without AS in the same gene
#only focused on SE and RI
#SE
AS_SE_gene_list <- list_for_ML_data_f(gene_having_AS_info_l_f$exon, c("RI", "SE")) 

bind_rows(AS_SE_gene_list) %>%
  select(AS) %>%
  distinct()
gene_having_AS_info_l_f$exon %>%
  group_by(anno) %>%
  mutate(feature_count = n()) %>%
  filter(feature_count > 2) %>%
  select(-feature_count) %>%
  ungroup() %>%
  left_join(all_AS_intersected_feature_raw, by = c("feature", "anno", "order")) %>%
  replace_na(list(trt = "no", AS = "no", PI = "no", event = "no")) %>%
  #remove events in the first or second exon/intron (promoter effects)
  filter(order != 1 & order != 2) %>%
  mutate(AS_type_label = if_else(AS == "SE", 1, 0)) %>%
  group_by(anno) %>%
  mutate(sum_AS_type_label = sum(AS_type_label)) %>%
  filter(sum_AS_type_label != 0) %>%
  View()

#RI
AS_RI_gene_list <- list_for_ML_data_f(gene_having_AS_info_l_f$intron, c("SE", "RI")) 

bind_rows(AS_RI_gene_list) %>%
  select(AS) %>%
  distinct()
#create the function for random choosing the same amount of exon/intron as that of AS events occurred in that gene 
#(internal control)

min(unique(c(3,3,4,4,5,5)))
make_AS_ML_data_f <- function(data, a){ 
  data_a <- data %>%
    ungroup()
  
  AS_mat <- data_a %>%
    filter(PI == a[1]) 
  
  if(nrow(AS_mat)==0){
    noAS_mat <- tibble()
  } else {
    noAS_mat <- data_a %>%
      filter(AS == "no") %>%
      filter(!order %in% unique(c(c(unique(AS_mat$order)-1), unique(AS_mat$order), unique(AS_mat$order)+1)))
  }
  
  
  noAS_mat_f <- head(noAS_mat[sample(1:nrow(noAS_mat)),], n = nrow(AS_mat)) 
  
  AS_noAS_pair_mat <- bind_rows(AS_mat, noAS_mat_f) 
  AS_noAS_pair_mat
}


length(exon_data_for_ML_list)

bind_rows(exon_data_for_ML_list) %>%
  group_by(anno) %>%
  mutate(max_exon = max(order)) %>%
  select(anno, max_exon) %>%
  distinct() %>%
  group_by(max_exon) %>%
  summarise(sum_max_exon = n()) %>%
  View()
  

exon_data_for_ML_list[[4]] %>%
  ungroup() %>%
  filter(PI == "PSI_40_60")


make_AS_ML_data_f(AS_SE_gene_list$Solyc01g109080.3.MC3, "PSI_5_20")

#make list for parallelism 
PI_type = unique(all_AS_intersected_feature_raw$PI)

#create the list with randommly selected internal ctrl 
AS_SE_for_ML_list <- vector("list", length = length(PI_type))
AS_RI_for_ML_list <- vector("list", length = length(PI_type))

for (i in seq_along(PI_type)) {
  AS_SE_for_ML_list[[i]] <- AS_SE_gene_list %>%
    map(. %>% make_AS_ML_data_f(., PI_type[[i]]))
  
  AS_RI_for_ML_list[[i]] <- AS_RI_gene_list %>%
    map(. %>% make_AS_ML_data_f(., PI_type[[i]]))
}

length(AS_SE_for_ML_list[[1]])

#create the function for producing segments of AS events and internal ctrls
unique_segments_AS_and_intctrl_f <- function(data, a){
  data %>%
    bind_rows() %>%
    drop_na() %>%
    mutate(PI = if_else(AS == "no", a[1], PI)) %>%
    select(chr, str, end, strand, trt, AS, PI) %>%
    distinct() %>%
    mutate(seg_size = end - str, order_t = str_c(PI, "_segment_", c(1:n())))
}

AS_SE_for_ML_list[[1]] %>%
  bind_rows()  %>%
  drop_na() %>%
  mutate(PI = if_else(AS == "no", PI_type[[1]], PI))%>%
  select(chr, str, end, strand, trt, AS, PI) %>%
  distinct() %>%
  mutate(seg_size = end - str, order_t = str_c("_segment_", c(1:n())))

#create unique segments of AS events and their internal controls
unique_seg_AS_SE_for_ML_list <- vector("list", length = length(PI_type))
unique_seg_AS_RI_for_ML_list <- vector("list", length = length(PI_type))

for (i in seq_along(PI_type)) {
  unique_seg_AS_SE_for_ML_list[[i]] <-  unique_segments_AS_and_intctrl_f(AS_SE_for_ML_list[[i]], PI_type[[i]])
  
  unique_seg_AS_RI_for_ML_list[[i]] <- unique_segments_AS_and_intctrl_f(AS_RI_for_ML_list[[i]], PI_type[[i]])
}

#create the function for the output of internal controls

producing_output_internal_ctrl_f <- function(data){
  
  a <- data %>%
    filter(AS != "no")
  
  b <- data %>%
    filter(AS == "no") %>%
    select(chr, str, end, trt, AS, PI) %>%
    arrange(chr, str)
  
  write_delim(b, str_c("./data/AS_events_HS0/AS_HS0_PSI_segment/AS_HS0_", unique(a$AS), "_", 
                          unique(a$PI), "_in_ctrl.bed"), col_names = FALSE, delim = "\t")
}

#create the output of internal ctrl segments
lapply(unique_seg_AS_SE_for_ML_list, producing_output_internal_ctrl_f)

lapply(unique_seg_AS_RI_for_ML_list, producing_output_internal_ctrl_f)

#create the function for producing 100-bp extentions of AS events and their internal ctrls (for ML)
producing_100bp_extentions_f <- function(data){
bind_rows(data, data) %>%
  ungroup() %>%
  mutate(side = rep(c("five", "three"), each = nrow(data))) %>% 
  mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                         if_else(strand == "+" & side == "three", end,
                                 if_else(strand == "-" & side == "five", end, str - 100)))) %>%
  mutate(end_n = if_else(strand == "+" & side == "five", str,
                         if_else(strand == "+" & side == "three", end + 100,
                                 if_else(strand == "-" & side == "five", end + 100, str)))) %>%
  select(chr, str = str_n, end = end_n, strand, trt, AS, PI, seg_side = side, seg_size, order_t) %>%
  arrange(chr, str)
}

#produce two lists of 100-bp extentions
#SE
AS_SE_100bp_extention_l <- 
unique_seg_AS_SE_for_ML_list %>%
  map(. %>% producing_100bp_extentions_f)

#RI
AS_RI_100bp_extention_l <- 
  unique_seg_AS_RI_for_ML_list %>%
  map(. %>% producing_100bp_extentions_f)

#produce function for the output file
producing_output_100bp_extentions_f <- function(data){
  
  a <- data %>%
    filter(AS != "no")
  
  write_delim(data, str_c("./data/AS_bed_for_ML_HS0/segment_for_ML_", unique(a$trt), "_", 
                          unique(a$PI), "_", unique(a$AS), "_", "target_in_ctrl.bed"), col_names = FALSE, delim = "\t")
}

#write the output files
AS_SE_100bp_extention_l %>%
  map(. %>% producing_output_100bp_extentions_f)

lapply(AS_RI_100bp_extention_l, producing_output_100bp_extentions_f)

AS_RI_100bp_extention_l[[1]]
length(AS_RI_100bp_extention_l)
