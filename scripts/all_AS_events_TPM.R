#For the occurrence of alternative splicing, I am wondering which factors matters,
#gene size, which exon or intron tends to have larger frequency to have AS,
#The following analysis is for identifying these factors,
#And it would be helpful for selecting the control dataset for ploting the signal of histone marks
#Since there are two kinds of annotation data in our lab, 
#I first used SollycM82_genes_v1.1.0.gff3 (without any modification of non-coding RNA) to extract exons and introns
#Then, I focus on the other one "Sl_M82_20220504_all.gff3" with the modification of non-coding RNA by Jeremie

library(tidyverse)

#M82 without modifying non-coding RNA
M82_all_exon <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12Chr_exon.bed", 
                           delim = "\t", col_names = c("chr", "str", "end", "info", "strand"))
M82_all_intron <- read_delim("./data/M82_annotation_data/SollycM82_genes_v1.1.0_12Chr_intron.bed", 
                             delim = "\t", col_names = c("chr", "str", "end", "info", "strand"))

#get the location of exon and gene in "info" to extract the gene/exon name 
M82_exon_exon_str_loc <- str_locate(M82_all_exon$info, "ID=exon:")

M82_exon_gene_str_loc <- str_locate(M82_all_exon$info, "Parent=mRNA:")

M82_exon_gene_end_loc <- str_locate(M82_all_exon$info, ";extra_copy")
?str_remove()
M82_all_exon %>%
  mutate(gene_str_loc = M82_exon_gene_str_loc[,2] + 1, gene_end_loc = M82_exon_gene_end_loc[,1] - 1) %>%
  mutate(exon_str_loc = M82_exon_exon_str_loc[,2] + 1, exon_end_loc = M82_exon_gene_str_loc[,1] - 2) %>%
  mutate(gene = str_sub(info, gene_str_loc, gene_end_loc),
         exon = str_sub(info, exon_str_loc, exon_end_loc)) %>%
  mutate(gene = str_remove(gene, ";=")) %>%
  select(-gene_str_loc, -gene_end_loc, -exon_str_loc, -exon_end_loc) %>%
  #exon_label_a means the order from the M82 annotation data
  mutate(exon_label_a = str_remove(exon, str_c(gene, "."))) %>%
  group_by(gene) %>%
  #exon_label_m means the order done manually
  mutate(exon_label_m = if_else(strand == "+", c(1:n()), c(n():1))) %>%
  filter(exon_label_m == exon_label_a) %>%
  filter(strand == "-") 

#M82 with modifying non-coding RNA
#import data
all_AS_events <- read_delim("./data/rMATS_out/All_events_all_comparisons.tsv", delim = "\t") %>%
  mutate(ID_modified = str_remove(ID, GeneID))

#split ID_modified to get 6 position of the necessary information for alternative splicing
AS_coordinate_all <- str_split(all_AS_events$ID_modified, "_", simplify = TRUE)

#for all AS
all_AS_events_bed_raw_v1 <- tibble(
  chr = all_AS_events$chr,
  GeneID = all_AS_events$GeneID,
  strand = all_AS_events$strand,
  pos_1 = as.double(AS_coordinate_all[,2]),
  pos_2 = as.double(AS_coordinate_all[,3]),
  pos_3 = as.double(AS_coordinate_all[,4]),
  pos_4 = as.double(AS_coordinate_all[,5]),
  pos_5 = as.double(AS_coordinate_all[,6]),
  pos_6 = as.double(AS_coordinate_all[,7]),
  AS_type = all_AS_events$AS_type,
  comp = all_AS_events$comp,
  FDR = all_AS_events$FDR,
  IncLevel1 = all_AS_events$IncLevel1,
  IncLevel2 = all_AS_events$IncLevel2,
  Incdf = all_AS_events$IncLevelDifference
) 

#calculate TPM for all of these events
Inc_all_sample_1 <- str_split(all_AS_events$IJC_SAMPLE_1, ",", simplify = TRUE)
Sk_all_sample_1 <- str_split(all_AS_events$SJC_SAMPLE_1, ",", simplify = TRUE)
Inc_all_sample_2 <- str_split(all_AS_events$IJC_SAMPLE_2, ",", simplify = TRUE)
Sk_all_sample_2 <- str_split(all_AS_events$SJC_SAMPLE_2, ",", simplify = TRUE)
IncLevel_all_sample1_l <- str_split(all_AS_events_bed_raw_v1$IncLevel1, ",", simplify = TRUE)
IncLevel_all_sample2_l <- str_split(all_AS_events_bed_raw_v1$IncLevel2, ",", simplify = TRUE)



all_AS_events_bed_TPM_raw <- all_AS_events_bed_raw_v1 %>%
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
  select(chr, str, end, GeneID, strand, AS_type, comp, FDR, IncLevel1, IncLevel2, Incdf, GeneID, pos_1, pos_2, pos_3, pos_4, pos_5, pos_6) %>%
  mutate(Inc_S1_r1 = as.double(Inc_all_sample_1[,1]), Inc_S1_r2 = as.double(Inc_all_sample_1[,2]),
         Sk_S1_r1 = as.double(Sk_all_sample_1[,1]), Sk_S1_r2 = as.double(Sk_all_sample_1[,2]),
         Inc_S2_r1 = as.double(Inc_all_sample_2[,1]), Inc_S2_r2 = as.double(Inc_all_sample_2[,2]),
         Sk_S2_r1 = as.double(Sk_all_sample_2[,1]), Sk_S2_r2 = as.double(Sk_all_sample_2[,2])) %>%
  mutate(S1 = str_sub(comp, 1, 3), S2 = str_sub(comp, 8, 10)) %>%
  mutate(Inc_length = if_else(AS_type == "RI", pos_2 - pos_1,
                              if_else(AS_type == "SE", pos_6 - pos_5 + pos_4 - pos_3 + pos_2 - pos_1, pos_6 - pos_5 + pos_2 - pos_1))) %>%
  mutate(Sk_length = pos_6 - pos_5 + pos_4 - pos_3) %>%
  mutate(Inc_length = Inc_length/10^3, Sk_length = Sk_length/10^3) %>%
  mutate(Inclvl1_r1 = as.double(IncLevel_all_sample1_l[,1]), Inclvl1_r2 = as.double(IncLevel_all_sample1_l[,2]),
         Inclvl2_r1 = as.double(IncLevel_all_sample2_l[,1]), Inclvl2_r2 = as.double(IncLevel_all_sample2_l[,2])) %>%
  replace_na(list(Inclvl1_r1 = 10^5, Inclvl1_r2 = 10^5, Inclvl2_r1 = 10^5, Inclvl2_r2 = 10^5)) %>%
  group_by(comp) %>%
  mutate(count_S1_r1 = sum(Inc_S1_r1/Inc_length + Sk_S1_r1/Sk_length)/10^6,
         count_S1_r2 = sum(Inc_S1_r2/Inc_length + Sk_S1_r2/Sk_length)/10^6,
         count_S2_r1 = sum(Inc_S2_r1/Inc_length + Sk_S2_r1/Sk_length)/10^6,
         count_S2_r2 = sum(Inc_S2_r2/Inc_length + Sk_S2_r2/Sk_length)/10^6) %>%
  ungroup() %>%
  mutate(TPM_Inc_S1_r1 = Inc_S1_r1/Inc_length/count_S1_r1,
         TPM_Inc_S1_r2 = Inc_S1_r2/Inc_length/count_S1_r2,
         TPM_Sk_S1_r1 = Sk_S1_r1/Sk_length/count_S1_r1,
         TPM_Sk_S1_r2 = Sk_S1_r2/Sk_length/count_S1_r2,
         TPM_Inc_S2_r1 = Inc_S2_r1/Inc_length/count_S2_r1,
         TPM_Inc_S2_r2 = Inc_S2_r2/Inc_length/count_S2_r2,
         TPM_Sk_S2_r1 = Sk_S2_r1/Sk_length/count_S2_r1,
         TPM_Sk_S2_r2 = Sk_S2_r2/Sk_length/count_S2_r2) %>%
  mutate(TPM_S1_r1 = TPM_Inc_S1_r1 + TPM_Sk_S1_r1,
         TPM_S1_r2 = TPM_Inc_S1_r2 + TPM_Sk_S1_r2,
         TPM_S2_r1 = TPM_Inc_S2_r1 + TPM_Sk_S2_r1,
         TPM_S2_r2 = TPM_Inc_S2_r2 + TPM_Sk_S2_r2)

quantile(all_AS_events_bed_TPM_raw$TPM_S2_r2, 0.05)

produce_invervals_TPM_threshold <- function(a){
all_AS_events_bed_TPM_raw %>%
  filter(TPM_S1_r1 > quantile(TPM_S1_r1, a[1]) &
         TPM_S1_r2 > quantile(TPM_S1_r2, a[1]) &
           TPM_S2_r1 > quantile(TPM_S2_r1, a[1]) &
           TPM_S2_r2 > quantile(TPM_S2_r2, a[1]))
}

all_AS_events_bed_TPM_q05 <- produce_invervals_TPM_threshold(0.05)

quantile(all_AS_events_bed_TPM_q05$Inclvl1_r1, seq(0, 1, 0.05))  

#add PI to all AS event before assigning sig and no-sig AS
all_AS_events_bed_TPM_q05_raw <- all_AS_events_bed_TPM_q05 %>%
  #filter(FDR > 0.05) %>%
  #PI means percentage of Inclusion (H: PI > 0.8, MH: 0.6 < PI <= 0.8, M: 0.4 < PI <= 0.6, ML: 0.2 < PI <= 0.4, L PI <= 0.2)
  mutate(PI = if_else(Inclvl1_r1 > 0.8 & Inclvl1_r2 > 0.8 & Inclvl2_r1 > 0.8 & Inclvl2_r2 > 0.8, "PI_H",
                      if_else(Inclvl1_r1 <= 0.2 & Inclvl1_r2 <= 0.2 & Inclvl2_r1 <= 0.2 & Inclvl2_r2 <= 0.2, "PI_L",
                              if_else(Inclvl1_r1 > 0.6 & Inclvl1_r1 <= 0.8 & Inclvl1_r2 > 0.6 & Inclvl1_r2 <= 0.8 &
                                        Inclvl2_r1 > 0.6 & Inclvl2_r1 <= 0.8 & Inclvl2_r2 > 0.6 & Inclvl2_r2 <= 0.8, "PI_MH",
                                      if_else(Inclvl1_r1 > 0.4 & Inclvl1_r1 <= 0.6 & Inclvl1_r2 > 0.4 & Inclvl1_r2 <= 0.6 &
                                                Inclvl2_r1 > 0.4 & Inclvl2_r1 <= 0.6 & Inclvl2_r2 > 0.4 & Inclvl2_r2 <= 0.6, "PI_M",
                                              if_else(Inclvl1_r1 > 0.2 & Inclvl1_r1 <= 0.4 & Inclvl1_r2 > 0.2 & Inclvl1_r2 <= 0.4 &
                                                        Inclvl2_r1 > 0.2 & Inclvl2_r1 <= 0.4 & Inclvl2_r2 > 0.2 & Inclvl2_r2 <= 0.4, "PI_ML", "others")))))) %>%
  filter(PI != "others")

all_DAS_events_bed_TPM_q05_raw <- all_AS_events_bed_TPM_q05 %>%
  filter(FDR < 0.05) %>%
  filter(FDR < 0.05) %>%
  mutate(PI = if_else(Inclvl1_r1 > 0.8 & Inclvl1_r2 > 0.8, "PI_H",
                      if_else(Inclvl1_r1 <= 0.2 & Inclvl1_r2 <= 0.2, "PI_L",
                              if_else(Inclvl1_r1 > 0.6 & Inclvl1_r1 <= 0.8 & Inclvl1_r2 > 0.6 & Inclvl1_r2 <= 0.8, "PI_MH",
                                      if_else(Inclvl1_r1 > 0.4 & Inclvl1_r1 <= 0.6 & Inclvl1_r2 > 0.4 & Inclvl1_r2 <= 0.6, "PI_M",
                                              if_else(Inclvl1_r1 > 0.2 & Inclvl1_r1 <= 0.4 & Inclvl1_r2 > 0.2 & Inclvl1_r2 <= 0.4, "PI_ML", "others")))))) %>%
  filter(PI != "others") 


#create the plot of non-sig and sig ASs for the number of their counts
#AS
all_AS_events_bed_TPM_q05_control_count_p <- all_AS_events_bed_TPM_q05_raw %>%
  filter(FDR > 0.05) %>%
  group_by(comp, AS_type, PI) %>%
  summarise(count = n()) %>%
  group_by(comp, AS_type) %>%
  mutate(count_comp_AS = sum(count)) %>%
  mutate(sd = sqrt(count)) %>%
  ggplot() + 
  geom_bar(aes(x = PI, y = count), stat = 'identity') +
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd, x = PI), width=.1, position = position_dodge(width = 0.9)) +
  #coord_cartesian(ylim = c(21000, 23000)) +
  facet_grid(comp ~ AS_type) +
  scale_y_continuous(breaks=seq(0, 25000, 5000))+
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))

ggsave("./analysis_output/all_AS_events_bed_TPM_q05_control_count_p.jpeg", all_AS_events_bed_TPM_q05_control_count_p, width = 500, height = 400, units = c("mm"), dpi = 320)


#DAS                                  
all_DAS_events_bed_TPM_q05_count_p <- all_DAS_events_bed_TPM_q05_raw %>%
  filter(FDR < 0.05) %>%
  group_by(comp, AS_type, PI) %>%
  summarise(count = n()) %>%
  group_by(comp, AS_type) %>%
  mutate(count_comp_AS = sum(count)) %>%
  mutate(sd = sqrt(count)) %>%
  ggplot() + 
  geom_bar(aes(x = PI, y = count), stat = 'identity') +
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd, x = PI), width=.1, position = position_dodge(width = 0.9)) +
  facet_grid(comp ~ AS_type) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(size = 18, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), axis.text.x = element_text(size = 16, face = "bold"))

ggsave("./analysis_output/all_DAS_events_bed_TPM_q05_count_p.jpeg", all_DAS_events_bed_TPM_q05_count_p, width = 500, height = 400, units = c("mm"), dpi = 320)


#generate the the list of control set using 5 categories of PI (percentage of Inclusion)
all_AS_events_bed_TPM_q05_control_l <- all_AS_events_bed_TPM_q05_raw %>%
  filter(FDR > 0.05) %>%
  #group_by(comp, AS_type, PI) %>%
  select(chr, str, end, comp, AS_type, PI) %>%
  arrange(chr, str, end) %>%
  mutate(comp = str_replace(comp, ".vs.", "_")) %>%
  split(.$comp) %>%
  map(. %>% split(.$AS_type)) %>%
  map(. %>% map(. %>% split(.$PI)))

#Also generate the list of control set without categorizing PI
all_AS_events_bed_TPM_q05_control_l_noPI <- all_AS_events_bed_TPM_q05_control_l %>%
  map(. %>% map(. %>% bind_rows))

#write control files (with PI groups) 
for (i in seq_along(all_AS_events_bed_TPM_q05_control_l)) {
  for (j in seq_along(all_AS_events_bed_TPM_q05_control_l[[1]])) {
    for (k in seq_along(all_AS_events_bed_TPM_q05_control_l[[1]][[1]])) {
      
      a <- str_c("data/AS_control_set/AS_control_TPM_q05_", names(all_AS_events_bed_TPM_q05_control_l)[[i]], "_", 
            names(all_AS_events_bed_TPM_q05_control_l[[i]])[[j]], "_", names(all_AS_events_bed_TPM_q05_control_l[[i]][[j]])[[k]], ".bed")
      
      write_delim(all_AS_events_bed_TPM_q05_control_l[[i]][[j]][[k]], a, delim = "\t", col_names = FALSE)
      
      
    }
  }
}

#write control files (without PI groups)
for (i in seq_along(all_AS_events_bed_TPM_q05_control_l_noPI)) {
  for (j in seq_along(all_AS_events_bed_TPM_q05_control_l_noPI[[1]])) {
    
      a <- str_c("data/AS_control_set/AS_control_TPM_q05_", names(all_AS_events_bed_TPM_q05_control_l_noPI)[[i]], "_", 
                 names(all_AS_events_bed_TPM_q05_control_l_noPI[[i]])[[j]], ".bed")
      
      write_delim(all_AS_events_bed_TPM_q05_control_l_noPI[[i]][[j]], a, delim = "\t", col_names = FALSE)
      
      
    }
  }


str_c("data/AS_control_set/AS_control_TPM_q05_", names(all_AS_events_bed_TPM_q05_control_l_noPI)[[i]], "_", 
      names(all_AS_events_bed_TPM_q05_control_l_noPI[[i]])[[j]], ".bed")

#generate the list of sig event using FDR below 0.05 
#PI added
all_AS_events_bed_TPM_q05_FDR05_l <- 
  all_DAS_events_bed_TPM_q05_raw %>%
  select(chr, str, end, comp, AS_type, PI) %>%
  mutate(comp = str_replace(comp, ".vs.", "_")) %>%
  arrange(chr, str, end) %>%
  group_by(comp, AS_type, PI) %>%
  split(.$comp) %>%
  map(. %>% split(.$AS_type)) %>%
  map(. %>% map(. %>% split(.$PI)))

#without PI
all_AS_events_bed_TPM_q05_FDR05_l_noPI <- all_AS_events_bed_TPM_q05_FDR05_l %>%
  map(. %>% map(. %>% bind_rows))

names(all_AS_events_bed_TPM_q05_FDR05_l[[1]])


#write sig event with PI
for (i in seq_along(all_AS_events_bed_TPM_q05_FDR05_l)) {
  for (j in seq_along(all_AS_events_bed_TPM_q05_FDR05_l[[1]])) {
    for (k in seq_along(all_AS_events_bed_TPM_q05_FDR05_l[[1]][[1]])) {
      
    
    
    a <- str_c("data/AS_control_set/AS_sig_FDR05_TPM_q05_", names(all_AS_events_bed_TPM_q05_FDR05_l)[[i]], "_", 
               names(all_AS_events_bed_TPM_q05_FDR05_l[[i]])[[j]], "_", names(all_AS_events_bed_TPM_q05_FDR05_l[[i]][[j]])[[k]], ".bed")
    
    write_delim(all_AS_events_bed_TPM_q05_FDR05_l[[i]][[j]][[k]], a, delim = "\t", col_names = FALSE)
    
    }
  }
}

#write sig event without PI
for (i in seq_along(all_AS_events_bed_TPM_q05_FDR05_l_noPI)) {
  for (j in seq_along(all_AS_events_bed_TPM_q05_FDR05_l_noPI[[1]])) {
    
      a <- str_c("data/AS_control_set/AS_sig_FDR05_TPM_q05_", names(all_AS_events_bed_TPM_q05_FDR05_l_noPI)[[i]], "_", 
                 names(all_AS_events_bed_TPM_q05_FDR05_l_noPI[[i]])[[j]], ".bed")
      
      write_delim(all_AS_events_bed_TPM_q05_FDR05_l_noPI[[i]][[j]], a, delim = "\t", col_names = FALSE)
    
    
  }
}

#modify bed files of annotation data (Jeremie) to add order number of exon or intron
M82_rMATs_anno_all_exon <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_exon.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name"))
M82_rMATs_anno_all_intron <- read_delim("./data/M82_annotation_data/M82_rMATs_anno_all_intron.bed", delim = "\t", col_names = c("chr", "str", "end", "strand", "feature", "source", "gene_name"))


M82_rMATs_anno_all_feature <- bind_rows(M82_rMATs_anno_all_intron, M82_rMATs_anno_all_exon) %>%
  arrange(chr, str, end) %>%
  group_by(gene_name, feature) %>%
  mutate(order_n = if_else(strand == "+", c(1:n()), c(n():1))) %>%
  split(.$feature)

write_delim(M82_rMATs_anno_all_feature[[1]], "./data/M82_annotation_data/M82_rMATs_anno_all_exon_order.bed", col_names = FALSE, delim = "\t")
write_delim(M82_rMATs_anno_all_feature[[2]], "./data/M82_annotation_data/M82_rMATs_anno_all_intron_order.bed", col_names = FALSE, delim = "\t")


###check RNA levels of different events######
#check TPM between ctrl and sig set (no PI)
AS_NS_TPM_q05_ctrl_check_l <- all_AS_events_bed_TPM_q05_control_raw %>%
  #group_by(comp, AS_type) %>%
  #summarise(TPM_S1_r1_mean = median(TPM_S1_r1),
            #TPM_S1_r2_mean = median(TPM_S1_r2),
            #TPM_S2_r1_mean = median(TPM_S2_r1),
            #TPM_S2_r2_mean = median(TPM_S2_r2), count = n()) %>%
  split(.$AS_type) %>%
  map(. %>% split(.$comp))

AS_sig_TPM_q05_ctrl_check_l <- all_AS_events_bed_TPM_q05 %>%
  filter(FDR < 0.05) %>%
  mutate(PI = if_else(Inclvl1_r1 > 0.8 & Inclvl1_r2 > 0.8, "PI_H",
                      if_else(Inclvl1_r1 <= 0.2 & Inclvl1_r2 <= 0.2, "PI_L",
                              if_else(Inclvl1_r1 > 0.6 & Inclvl1_r1 <= 0.8 & Inclvl1_r2 > 0.6 & Inclvl1_r2 <= 0.8, "PI_MH",
                                      if_else(Inclvl1_r1 > 0.4 & Inclvl1_r1 <= 0.6 & Inclvl1_r2 > 0.4 & Inclvl1_r2 <= 0.6, "PI_M",
                                              if_else(Inclvl1_r1 > 0.2 & Inclvl1_r1 <= 0.4 & Inclvl1_r2 > 0.2 & Inclvl1_r2 <= 0.4, "PI_ML", "others")))))) %>%
  filter(PI != "others") %>%
  #group_by(comp, AS_type, PI) %>%
  #summarise(TPM_S1_r1_mean = median(TPM_S1_r1),
            #TPM_S1_r2_mean = median(TPM_S1_r2),
            #TPM_S2_r1_mean = median(TPM_S2_r1),
            #TPM_S2_r2_mean = median(TPM_S2_r2),
            #count = n()) %>%
  split(.$AS_type) %>%
  map(. %>% split(.$comp))

tibble(TPM_S1_r1_NS = )

t.test(AS_NS_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S1_r1, AS_sig_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S1_r2)
t.test(AS_NS_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S2_r1, AS_NS_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S2_r2)


?aov
mean(AS_NS_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S2_r2)
mean(AS_sig_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0$TPM_S2_r2)
?spread
AS_NS_RI_6_0_check <- AS_NS_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0 %>%
  select(TPM_S1_r1, TPM_S1_r2, TPM_S2_r1, TPM_S2_r2) %>%
  mutate(AS_status = "NS") %>%
  gather(key = "data", value = "TPM_value", c(1:4)) %>%
  mutate(sample = str_c(data, "_", AS_status)) %>%
  select(sample, TPM_value)

AS_sig_RI_6_0_check <- AS_sig_TPM_q05_ctrl_check_l$RI$HS6.vs.HS0 %>%
  select(TPM_S1_r1, TPM_S1_r2, TPM_S2_r1, TPM_S2_r2) %>%
  mutate(AS_status = "sig") %>%
  gather(key = "data", value = "TPM_value", c(1:4)) %>%
  mutate(sample = str_c(data, "_", AS_status)) %>%
  select(sample, TPM_value)

check_RNA_level_RI_6_0_raw <- bind_rows(AS_NS_RI_6_0_check, AS_sig_RI_6_0_check)

check_RNA_level_RI_6_0_raw_S1 <- check_RNA_level_RI_6_0_raw %>%
  filter(str_detect(sample, "S1")==TRUE)

check_RNA_level_RI_6_0_raw_S2 <- check_RNA_level_RI_6_0_raw %>%
  filter(str_detect(sample, "S2")==TRUE)

aov_RNA_level_RI_6_0_raw_S1 <- aov(TPM_value ~ sample, data = check_RNA_level_RI_6_0_raw_S1)
aov_RNA_level_RI_6_0_raw_S2 <- aov(TPM_value ~ sample, data = check_RNA_level_RI_6_0_raw_S2)

summary(aov_RNA_level_RI_6_0_raw_S2)
