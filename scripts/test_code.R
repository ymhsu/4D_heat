all_AS_intersected_feature_raw %>% 
  left_join(M82_rMATs_anno_all_gene_light) %>%
  mutate(gene_size = gene_size %/% 1000) %>%
  group_by(gene_size, AS, event, feature, order) %>%
  summarise(count = n()) %>%
  View()

all_AS_intersected_feature_raw %>%
  group_by(AS, event, feature, order) %>%
  summarise(count = n()) %>%
  View()

all_AS_intersected_feature_raw %>%
  filter(anno == "Solyc_chr1_NPC00337.MC1")
gene_having_AS_info_l_f$exon %>%
  filter(anno == "Solyc_chr1_NPC00337.MC1")
M82_rMATs_anno_all_exon_order %>%
  group_by(anno) %>%
  mutate(pre_end = lag(end)) %>%
  drop_na() %>%
  mutate(dist_exon = str - pre_end) %>%
  ungroup() %>%
  summarise(q25 = quantile(dist_exon, 0.25), q50 = median(dist_exon), q75 = quantile(dist_exon, 0.75))


exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2

make_AS_ML_data_f(exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2, "HS6_HS0", "PI_H", "RI")

exon_intron_data_for_ML_list$exon$Solyc_chr2_NPC00216.MC2 %>%
  filter(event == "no") %>%
  filter(order != 4-1 & order != 4+1)

make_AS_ML_data_f <- function(data, a, b, c){ 
data_a <- data %>%
  ungroup()
AS_sig_m <- data_a %>%
  filter(comp == a[1] & PI == b[1] &  AS == c[1] & event == "sig") 

AS_non_sig_m <- data_a %>%
  filter(comp == a[1] & PI == b[1] &  AS == c[1] & event == "non-sig") 

if(nrow(AS_sig_m)==0){
  AS_noAS_sig_m <- tibble()
} else {
  AS_noAS_sig_m <- data_a %>%
    filter(event == "no") %>%
    filter(!order %in% c(c(unique(AS_sig_m$order)-1), unique(AS_sig_m$order), unique(AS_sig_m$order)+1))
}

if(nrow(AS_non_sig_m)==0){
  AS_noAS_non_sig_m <- tibble()
} else {
  AS_noAS_non_sig_m <- data_a %>%
    filter(event == "no") %>%
    filter(!order %in% c(c(unique(AS_non_sig_m$order)-1), unique(AS_non_sig_m$order), unique(AS_non_sig_m$order)+1))
}

AS_noAS_sig_m_f <- head(AS_noAS_sig_m[sample(1:nrow(AS_noAS_sig_m)),], n = nrow(AS_sig_m)) 

AS_sig_pair_f <- bind_rows(AS_sig_m, AS_noAS_sig_m_f) %>%
  mutate(pair = "sig")

AS_noAS_non_sig_m_f <- head(AS_noAS_non_sig_m[sample(1:nrow(AS_noAS_non_sig_m)),], n = nrow(AS_non_sig_m)) 

AS_non_sig_pair_f <- bind_rows(AS_non_sig_m, AS_noAS_non_sig_m_f) %>%
  mutate(pair = "non-sig")

bind_rows(AS_sig_pair_f, AS_non_sig_pair_f)
}
length(exon_intron_data_for_ML_list$exon)
ML_data_test_head10 <- list(head(exon_intron_data_for_ML_list$exon, n = 10), head(exon_intron_data_for_ML_list$intron, n = 10))
system.time(ML_data_test_head10 %>%
  map(. %>% map(. %>% make_AS_ML_data_f(., comb_comp_PI_AS_list[[1]][[1]]))) %>%
  map(. %>% bind_rows()))

exon_intron_data_for_ML_comb_l[[1]][[1]][[4]][[4]] <- ML_data_test_head10[[1]] %>%
  map(. %>% make_AS_ML_data_f(., comp[[1]], PI_type[[4]], AS_type[[4]])) %>%
  bind_rows()

exon_intron_data_for_ML_list

registerDoParallel(cores = 10)
getDoParWorkers()


exon_intron_data_for_ML_list_v4$HS1_HS0$PI_H$A3S
bind_rows(exon_intron_data_for_ML_list_v4$HS1_HS0$PI_H$A3S, exon_intron_data_for_ML_list_v4$HS1_HS0$PI_H$A3S) %>%
  ungroup() %>%
  mutate(side = rep(c("five", "three"), each = nrow(exon_intron_data_for_ML_list_v4$HS1_HS0$PI_H$A3S))) %>% 
  mutate(str_n = if_else(strand == "+" & side == "five", str - 100,
                         if_else(strand == "+" & side == "three", end,
                                 if_else(strand == "-" & side == "five", end, str - 100)))) %>%
  mutate(end_n = if_else(strand == "+" & side == "five", str,
                         if_else(strand == "+" & side == "three", end + 100,
                                 if_else(strand == "-" & side == "five", end + 100, str)))) %>%
  select(chr, str = str_n, end = end_n, strand, feature, source, anno, order, comp, AS, PI, event, seg_side = side) %>%
  arrange(chr, str)



#produce data for ML
#make a list for operating 10 cores
#I will use 10 cores, so here I made the length of list as 12 (120/10 = 12)

AS_exon_intron_sig_no_event_l_10cores <- vector("list", length = length(AS_exon_intron_sig_no_event_l)/10)
names_AS_exon_intron_sig_no_event_l_10cores <- vector("list", length = length(AS_exon_intron_sig_no_event_l)/10)


for (i in seq_along(AS_exon_intron_sig_no_event_l_10cores)) {
  print(seq(10*(i-1)+1, 10*(i-1)+10, 1))
  #comb_comp_PI_AS_list[[i]] <- comb_comp_PI_AS_list_raw[seq(12*(i-1)+1, 12*(i-1)+12, 1)]
  AS_exon_intron_sig_no_event_l_10cores[[i]] <- AS_exon_intron_sig_no_event_l[seq(10*(i-1)+1, 10*(i-1)+10, 1)]
  names_AS_exon_intron_sig_no_event_l_10cores[[i]] <- names_AS_exon_intron_sig_no_event_l[seq(10*(i-1)+1, 10*(i-1)+10, 1)]
}

names_AS_exon_intron_sig_no_event_l_10cores[[1]][[1]]

registerDoParallel(cores = 10)
getDoParWorkers()

AS_all_seg_his_signal_raw <-
  foreach(i=c(1:19), .packages = c("tidyverse")) %:%
  foreach(j=c(1:12), .packages = c("tidyverse")) %:%
  foreach(k=1:10, .packages = c("tidyverse")) %dopar% {
    
    feature_size <- AS_exon_intron_sig_no_event_l_10cores[[j]][[k]] %>%
      group_by(anno, feature, order) %>%
      mutate(str_n = min(end), end_n = max(str)) %>%
      mutate(feature_size = end_n - str_n) %>%
      ungroup() %>%
      select(anno, feature, order, feature_size) %>%
      distinct()
    
    test <-
      read_delim(
        str_c(
          "./data/AS_bed_for_ML/all_segments/",
          histone_mark_list_m$name_raw[[i]],
          "_",
          str_remove(names_AS_exon_intron_sig_no_event_l_10cores[[j]][[k]], "segment_"),
          "_12chr.bed.gz"
        ),
        delim = "\t",
        col_names = c(
          "chr",
          "str",
          "end",
          "strand",
          "feature",
          "source",
          "anno",
          "order",
          "comp",
          "AS",
          "PI",
          "event",
          "pair",
          "seg_side",
          "AS_signal"
        )
      ) %>%
      group_by(anno, feature, order, comp, AS, PI, event, seg_side) %>%
      summarise(mean_signal = sum(AS_signal) / n()) %>%
      mutate(read_count = histone_mark_list_m$read_count[[i]]) %>%
      mutate(RPM = mean_signal / read_count * 10 ^ 6) %>%
      left_join(feature_size) %>%
      mutate(RPKM = RPM / feature_size) %>%
      select(-mean_signal, -read_count, -RPM, -feature_size)
    
    names(test) <-
      c(
        "anno",
        "order",
        "comp",
        "AS",
        "PI",
        "event",
        "seg_side",
        "feature",
        str_c(histone_mark_list_m$name_light[[i]], "_RPKM")
      )
    
    AS_exon_intron_sig_no_event_l_10cores[[j]][[k]] %>%
      left_join(test)
  }

AS_all_seg_his_signal_raw_v2 <- vector("list", length = length(AS_exon_intron_sig_no_event_l))

for (i in seq_along(histone_mark_list_m$name_raw)) {
  for (j in seq_along(AS_exon_intron_sig_no_event_l_10cores)) {
    for (k in seq_along(AS_exon_intron_sig_no_event_l_10cores[[1]])) {
      print(c(10*(j-1)+k))
      AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]] <- append(AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]], list(AS_all_seg_his_signal_raw[[i]][[j]][[k]]))
      AS_all_seg_his_signal_raw_v2[[c(10*(j-1)+k)]] %>%
        reduce(full_join)
    }
  }
}    


AS_all_seg_his_signal_raw_v2 <- AS_all_seg_his_signal_raw_v2 %>% 
  map(. %>% reduce(full_join))

output_name_f <- function(data){
  
  if(data$event[[1]]=="no") {
    output_name <- str_c("table_for_ML_ready_", data$comp[[1]], "_", data$PI[[1]], "_", data$AS[[1]], "_", "ctrl")
  } else {
    output_name <- str_c("table_for_ML_ready_", data$comp[[1]], "_", data$PI[[1]], "_", data$AS[[1]], "_", "target")
  }  
  
  output_name
}


for (i in seq_along(AS_all_seg_his_signal_raw_v2)) {
  out_name <- output_name_f(AS_all_seg_his_signal_raw_v2[[i]])
  
  write_delim(AS_all_seg_his_signal_raw_v2[[i]], str_c("./data/AS_bed_for_ML/all_segments/", out_name), delim = "\t", col_names = TRUE)
}
version

output_name_f(AS_all_seg_his_signal_raw_v2[[2]])

install.packages("Boruta")
library(Boruta)


table_ready_for_ML_raw_ctrl <- table_ready_for_ML_raw[which(str_detect(name_table_ready_for_ML_raw, "ctrl"))]
names(table_ready_for_ML_raw_ctrl)

test_H1_H0_H_L_A5S <- bind_rows(table_ready_for_ML_raw_ctrl[[21]], table_ready_for_ML_raw_ctrl[[29]]) %>%
  select(-comp, -AS, -PI) %>%
  as.data.frame()

test_H1_H0_L_M_A5S <- bind_rows(table_ready_for_ML_raw_ctrl[[37]], table_ready_for_ML_raw_ctrl[[29]]) %>%
  select(-comp, -AS, -PI) %>%
  as.data.frame()

test_H1_H0_H_L_A5S_ctrl_boruta <- Boruta(type~., test_H1_H0_H_L_A5S, doTrace = 2)
attStats(test_H1_H0_H_L_A5S_ctrl_boruta)

test_H1_H0_L_M_A5S_ctrl_boruta <- Boruta(type~., test_H1_H0_L_M_A5S, doTrace = 2)

process_the_Boruta_data(test_H1_H0_H_L_A5S_ctrl_boruta) %>%
  pivot_longer(everything()) %>%
  mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
  ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
  geom_boxplot(aes(fill = mark)) +
  #ggtitle(name_ML_stage_1_target_for_plot[[i]]) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5))

process_the_Boruta_data(test_H1_H0_L_M_A5S_ctrl_boruta) %>%
  pivot_longer(everything()) %>%
  mutate(mark = if_else(str_detect(name, "shadow")==TRUE, "shadowmarks", "marks")) %>%
  ggplot(aes(x = fct_reorder(name, value, median), y = value)) +
  geom_boxplot(aes(fill = mark)) +
  #ggtitle(name_ML_stage_1_target_for_plot[[i]]) +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5))




#test_code
#for (i in seq_along(AS_exon_intron_sig_no_event_l_peak)) {
  #for (j in seq_along(histone_mark_list_m$name_raw)) {

#use RPKM to normalize the signal of 19 features for all samples
#RPKM = numReads / ( geneLength/1000 * totalNumReads/1,000,000 ) (from: https://www.metagenomics.wiki/pdf/qc/RPKM)
#geneLength is replaced by the exon or intron length

#create the function for the list of length of each intron/exon
feature_size_f <- function(data){
   data %>%
    group_by(anno, feature, order) %>%
    mutate(str_n = min(end), end_n = max(str)) %>%
    mutate(feature_size = end_n - str_n) %>%
    ungroup() %>%
    select(anno, feature, order, seg_side, order_t, feature_size) %>%
    distinct() %>%
    mutate(order_t = c(1:n()))
}

#make the list of feature size for all target or ctrl segments
AS_exon_intron_sig_no_event_size_l <- 
AS_exon_intron_sig_no_event_l %>% 
  map(., feature_size_f)
test <- 
bind_rows(table_ready_for_ML_raw_peak[[1]], table_ready_for_ML_raw_peak[[2]]) %>%
  replace_na(list(Pol2_RPKM = 0, H3K14ac_RPKM = 0, H3K4ac_RPKM = 0, K18ac_RPKM = 0, K27ac_RPKM = 0,
                  K27me3_RPKM = 0, K4me1_RPKM = 0, K4me3_RPKM = 0, K9ac_RPKM = 0, H3K36me3_RPKM = 0,
                  H3K79ac_RPKM = 0, H4K12ac_RPKM = 0, H4K16ac_RPKM = 0, H4K20ac_RPKM = 0, H4K5ac_RPKM = 0,
                  H4K8ac_RPKM = 0, K36ac_RPKM = 0, K4me2_RPKM = 0, K56ac_RPKM = 0)) %>%
select(-comp, -AS, -PI) %>%
  as.data.frame()

test_boruta <- Boruta(type~., test, doTrace = 2)
plot(test_boruta)



##plots for importance feature investigation
important_feature_f <- function(data, data2){
  process_the_Boruta_data(data) %>%
    as_tibble() %>%
    gather(key = feature, value = feature_value) %>%
    group_by(feature) %>%
    summarise(mean_feature_value = mean(feature_value)) %>%
    arrange(-mean_feature_value) %>%
    head(., 5) %>%
    mutate(name = data2)
}


ML_stage_1_peak_im_feature_raw <- vector("list", length = length(name_ML_stage_1_target_peaks_for_plot))


for (i in seq_along(name_ML_stage_1_target_for_plot)) {
  ML_stage_1_peak_im_feature_raw[[i]] <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1_peaks/", name_ML_stage_1_target_peaks_for_plot[[i]], "_peaks"), delim = "\t", col_names = TRUE) %>%
    gather(key = feature, value = feature_value) %>%
    group_by(feature) %>%
    summarise(mean_feature_value = mean(feature_value)) %>%
    arrange(-mean_feature_value) %>%
    head(., 5) %>%
    mutate(name = name_ML_stage_1_target_peaks_for_plot[[i]])
  
}

ML_stage_1_peak_im_feature_list <- bind_rows(ML_stage_1_peak_im_feature_raw) %>%
  mutate(comp = if_else(str_detect(name, "HS1_HS0")==TRUE, "HS1_HS0",
                        if_else(str_detect(name, "HS6_HS0")==TRUE, "HS6_HS0", "HS1_HS6"))) %>%
  mutate(AS = if_else(str_detect(name, "A5S")==TRUE, "A5S",
                      if_else(str_detect(name, "A3S")==TRUE, "A3S",
                              if_else(str_detect(name, "SE")==TRUE, "SE", "RI")))) %>%
  mutate(PI_comp = str_remove(name, str_c(comp, "_"))) %>%
  mutate(PI_comp = str_remove(PI_comp, str_c("_", AS, "_target"))) %>%
  select(-name)

ML_stage_1_peak_im_feature_list %>%
  group_by(comp, AS, feature) %>%
  summarise(count = n()) %>%
  arrange(comp, AS, -count) %>%
  filter(comp == "HS1_HS0") %>%
  ggplot() +
  geom_bar(aes(feature, count), stat = "identity") +
  facet_wrap(~AS, scales = "free") +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5))

#check the important features of boruta result from all signals
ML_stage_1_im_feature_raw <- vector("list", length = length(name_ML_stage_1_target_for_plot))


for (i in seq_along(name_ML_stage_1_target_for_plot)) {
  ML_stage_1_im_feature_raw[[i]] <- read_delim(str_c("./data/AS_bed_for_ML/all_segments/ML_results_table_stage_1/", name_ML_stage_1_target_for_plot[[i]]), delim = "\t", col_names = TRUE) %>%
    gather(key = feature, value = feature_value) %>%
    group_by(feature) %>%
    summarise(mean_feature_value = mean(feature_value)) %>%
    arrange(-mean_feature_value) %>%
    head(., 5) %>%
    mutate(name = name_ML_stage_1_target_for_plot[[i]])
  
}

ML_stage_1_im_feature_list <- bind_rows(ML_stage_1_im_feature_raw) %>%
  mutate(comp = if_else(str_detect(name, "HS1_HS0")==TRUE, "HS1_HS0",
                        if_else(str_detect(name, "HS6_HS0")==TRUE, "HS6_HS0", "HS1_HS6"))) %>%
  mutate(AS = if_else(str_detect(name, "A5S")==TRUE, "A5S",
                      if_else(str_detect(name, "A3S")==TRUE, "A3S",
                              if_else(str_detect(name, "SE")==TRUE, "SE", "RI")))) %>%
  mutate(PI_comp = str_remove(name, str_c(comp, "_"))) %>%
  mutate(PI_comp = str_remove(PI_comp, str_c("_", AS, "_target"))) %>%
  select(-name)

ML_stage_1_im_feature_list %>%
  group_by(comp, AS, feature) %>%
  summarise(count = n()) %>%
  arrange(comp, AS, -count) %>%
  filter(comp == "HS1_HS0") %>%
  ggplot() +
  geom_bar(aes(feature, count), stat = "identity") +
  facet_wrap(~AS, scales = "free") +
  theme(strip.text.x = element_text(colour = "black", face = "bold", size = 18), legend.text = element_text(size = 12, face = "bold"), plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 18),
        legend.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), 
        axis.text.y = element_text(color = "black", size = 12, face = "bold"), strip.text.y = element_text(colour = "black", face = "bold", size = 18), 
        axis.text.x = element_text(colour = "black", size = 10, face = "bold", angle = 45, vjust = 0.5))

attStats(table_ready_for_ML_stage_1_peak_ctrl_boruta[[10]]) %>%
  arrange(-meanImp)

name_pairwise_table_for_ML_stage_1_raw_peak_ctrl[[10]]



#make_AS_ML_outside_ctrl_data_f <- function(data, a){ 
  data_a <- data %>%
    ungroup()
  AS_sig_m <- data_a %>%
    filter(comp == a[1] & PI == a[2] &  AS == a[3] & event == "sig") 
  
  AS_non_sig_m <- data_a %>%
    filter(comp == a[1] & PI == a[2] &  AS == a[3] & event == "non-sig") 
  
  if(nrow(AS_sig_m)==0){
    AS_noAS_sig_m <- tibble()
  } else {
    AS_noAS_sig_m <- data_a %>%
      filter(event == "no") %>%
      filter(!order %in% c(c(unique(AS_sig_m$order)-1), unique(AS_sig_m$order), unique(AS_sig_m$order)+1))
  }
  
  if(nrow(AS_non_sig_m)==0){
    AS_noAS_non_sig_m <- tibble()
  } else {
    AS_noAS_non_sig_m <- data_a %>%
      filter(event == "no") %>%
      filter(!order %in% c(c(unique(AS_non_sig_m$order)-1), unique(AS_non_sig_m$order), unique(AS_non_sig_m$order)+1))
  }
  
  AS_noAS_sig_m_f <- head(AS_noAS_sig_m[sample(1:nrow(AS_noAS_sig_m)),], n = nrow(AS_sig_m)) 
  
  AS_sig_pair_f <- bind_rows(AS_sig_m, AS_noAS_sig_m_f) %>%
    mutate(pair = "sig")
  
  AS_noAS_non_sig_m_f <- head(AS_noAS_non_sig_m[sample(1:nrow(AS_noAS_non_sig_m)),], n = nrow(AS_non_sig_m)) 
  
  AS_non_sig_pair_f <- bind_rows(AS_non_sig_m, AS_noAS_non_sig_m_f) %>%
    mutate(pair = "non-sig")
  
  bind_rows(AS_sig_pair_f, AS_non_sig_pair_f)
#}

exon_intron_data_for_ML_list$exon$Solyc_chr1_NPC00216.MC1
gene_having_AS_info_l$exon
all_AS_intersected_feature_raw %>%
  group_by(comp, AS, PI) %>%
  summarise(count = n())
exon_intron_data_for_ML_list_v4$HS1_HS0$PI_H$A3S %>%
  left_join(M82_rMATs_anno_all_gene_light) %>%
  mutate(gene_size_dis_500bp = gene_size %/% 500) %>%
  select(anno, gene_size_dis_500bp) %>%
  distinct() %>%
  group_by(gene_size_dis_500bp) %>%
  summarise(count = n()) %>%
  View()
  
