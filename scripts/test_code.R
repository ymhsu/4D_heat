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
